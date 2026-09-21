"""Open Babel SDS-style PDBs must not gain a sulfate hydrogen."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

_TESTS_DIR = Path(__file__).resolve().parent
if str(_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(_TESTS_DIR))
import _paths  # noqa: E402


def setUpModule():
    _paths.ensure_ligandparam()


# Anionic sulfate oxygen as Open Babel writes it (element "O1-") plus a
# duplicate CONECT for S=O. 4 atoms: C-O-S-O-
_SDS_TAIL_PDB = """\
HETATM    1  C   UNL     1       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  O   UNL     1       1.400   0.000   0.000  1.00  0.00           O
HETATM    3  S   UNL     1       2.100   0.000   1.200  1.00  0.00           S
HETATM    4  O   UNL     1       3.500   0.000   1.200  1.00  0.00           O1-
CONECT    1    2
CONECT    2    1    3
CONECT    3    2    4    4
CONECT    4    3    3
END
"""


class TestPdbSanitize(unittest.TestCase):
    def test_o1minus_becomes_oxygen_with_charge(self):
        from ligandparam.io.Coordinates import _rewrite_pdb_atom_line

        line = (
            "HETATM   17  O   UNL     1      13.945   5.161   6.844  1.00  0.00"
            "           O1-\n"
        )
        out = _rewrite_pdb_atom_line(line)
        self.assertEqual(out[76:78], " O")
        self.assertEqual(out[78:80], "1-")
        self.assertEqual(len(out.rstrip("\n")), 80)

    def test_sanitize_keeps_atom_count_and_unique_conect(self):
        from ligandparam.io.Coordinates import count_structure_atoms, sanitize_pdb_ligand

        with tempfile.TemporaryDirectory() as td:
            src = Path(td) / "sds.pdb"
            dst = Path(td) / "sds.sanitized.pdb"
            src.write_text(_SDS_TAIL_PDB, encoding="utf-8")
            sanitize_pdb_ligand(src, dst)
            self.assertEqual(count_structure_atoms(src), 4)
            self.assertEqual(count_structure_atoms(dst), 4)
            text = dst.read_text(encoding="utf-8")
            o_line = [ln for ln in text.splitlines() if ln.startswith("HETATM") and ln.split()[1] == "4"][0]
            self.assertEqual(o_line[76:78], " O")
            self.assertEqual(o_line[78:80], "1-")
            conect = [ln for ln in text.splitlines() if ln.startswith("CONECT")]
            s_line = [ln for ln in conect if ln.split()[1] == "3"][0]
            partners = s_line.split()[2:]
            self.assertEqual(partners, ["2", "4"])

    def test_split_element_charge_tokens(self):
        from ligandparam.io.Coordinates import _split_pdb_element_charge

        self.assertEqual(_split_pdb_element_charge("O1-"), ("O", "1-"))
        self.assertEqual(_split_pdb_element_charge("O-1"), ("O", "1-"))
        self.assertEqual(_split_pdb_element_charge("N+"), ("N", "1+"))
        self.assertEqual(_split_pdb_element_charge("Cl"), ("Cl", ""))
        self.assertEqual(_split_pdb_element_charge("C"), ("C", ""))

    def test_drop_extra_hydrogen_against_reference(self):
        import numpy as np

        from ligandparam.io.Coordinates import match_current_to_reference

        ref = [
            ("C", np.array([0.0, 0.0, 0.0])),
            ("O", np.array([1.4, 0.0, 0.0])),
            ("S", np.array([2.1, 0.0, 1.2])),
            ("O", np.array([3.5, 0.0, 1.2])),
            ("H", np.array([-0.9, 0.0, 0.0])),
        ]
        extra_h = ("H", np.array([4.1, 0.4, 0.8]))
        current = list(ref) + [extra_h]
        elems, coords = match_current_to_reference(ref, current)
        self.assertEqual(len(elems), 5)
        self.assertEqual(elems.count("H"), 1)
        self.assertEqual(len(coords), 5)


# Methyl sulfate anion CH3OSO3- (9 atoms) plus a spurious sulfate H (10th).
_METHYL_SULFATE_MOL2 = """\
@<TRIPOS>MOLECULE
MSU
 10 9 1 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 C1           0.0000    0.0000    0.0000 C.3     1 MSU       0.0000
      2 H1          -0.5000    0.9000    0.0000 H       1 MSU       0.0000
      3 H2          -0.5000   -0.4500    0.7800 H       1 MSU       0.0000
      4 H3          -0.5000   -0.4500   -0.7800 H       1 MSU       0.0000
      5 OS           1.4300    0.0000    0.0000 O.3     1 MSU       0.0000
      6 S            2.9000    0.0000    0.0000 S.o2    1 MSU       0.0000
      7 O1           3.9000    1.1000    0.0000 O.3     1 MSU       0.0000
      8 HO           4.5000    1.5500    0.0000 H       1 MSU       0.0000
      9 O2           3.9000   -0.5500    0.9500 O.2     1 MSU       0.0000
     10 O3           3.9000   -0.5500   -0.9500 O.2     1 MSU       0.0000
@<TRIPOS>BOND
     1    1    2 1
     2    1    3 1
     3    1    4 1
     4    1    5 1
     5    5    6 1
     6    6    7 1
     7    7    8 1
     8    6    9 2
     9    6   10 2
"""


class TestMol2AnionSanitize(unittest.TestCase):
    def test_drops_sulfate_oh_when_charge_is_minus_one(self):
        from ligandparam.io.Coordinates import (
            closed_shell_ok,
            count_structure_atoms,
            parse_structure_atoms,
            sanitize_mol2_ligand,
        )

        with tempfile.TemporaryDirectory() as td:
            src = Path(td) / "msu.mol2"
            dst = Path(td) / "msu.sanitized.mol2"
            src.write_text(_METHYL_SULFATE_MOL2, encoding="utf-8")
            self.assertEqual(count_structure_atoms(src), 10)
            sanitize_mol2_ligand(src, dst, net_charge=-1)
            self.assertEqual(count_structure_atoms(dst), 9)
            elems = [e for e, _ in parse_structure_atoms(dst)]
            self.assertEqual(elems.count("H"), 3)
            self.assertTrue(closed_shell_ok(elems, charge=-1, multiplicity=1))
            text = dst.read_text(encoding="utf-8")
            self.assertIn("O.co2", text)

    def test_keeps_sulfate_oh_when_neutral(self):
        from ligandparam.io.Coordinates import count_structure_atoms, sanitize_mol2_ligand

        with tempfile.TemporaryDirectory() as td:
            src = Path(td) / "msu.mol2"
            dst = Path(td) / "msu.sanitized.mol2"
            src.write_text(_METHYL_SULFATE_MOL2, encoding="utf-8")
            sanitize_mol2_ligand(src, dst, net_charge=0)
            self.assertEqual(count_structure_atoms(dst), 10)

    def test_gaussian_atoms_drop_even_if_initial_mol2_is_protonated(self):
        from ligandparam.io.Coordinates import atoms_for_gaussian, count_structure_atoms

        with tempfile.TemporaryDirectory() as td:
            cwd = Path(td)
            protonated = cwd / "SDS.centered.mol2"
            protonated.write_text(_METHYL_SULFATE_MOL2, encoding="utf-8")
            (cwd / "SDS.initial.mol2").write_text(_METHYL_SULFATE_MOL2, encoding="utf-8")
            elems, coords = atoms_for_gaussian(
                protonated, cwd=cwd, net_charge=-1, multiplicity=1
            )
            self.assertEqual(len(elems), 9)
            self.assertEqual(len(coords), 9)
            self.assertEqual(count_structure_atoms(protonated), 10)
