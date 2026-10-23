import shutil
import tempfile
import unittest
from pathlib import Path

from ligandparam.stages.leap import StageLeap, _read_atom_type_elements

# 1-bromo-2-chloroethane, GAFF2 types from antechamber
MOL2 = """@<TRIPOS>MOLECULE
LIG
    8     7     1     0     0
SMALL
gas


@<TRIPOS>ATOM
      1 CL1         -1.9000    -0.9130    -0.6440 cl         1 LIG      -0.125820
      2 C1          -0.7140     0.2280    -0.0420 c3         1 LIG       0.032013
      3 C2           0.7040    -0.3200    -0.0560 c3         1 LIG       0.016740
      4 BR1          1.8130     1.1270     0.6530 br         1 LIG      -0.091526
      5 H1          -0.7930     1.1740    -0.5870 h1         1 LIG       0.044399
      6 H2          -0.9440     0.4670     1.0320 h1         1 LIG       0.044399
      7 H3           1.0420    -0.6270    -1.0460 h1         1 LIG       0.039898
      8 H4           0.7930    -1.1360     0.6890 h1         1 LIG       0.039898
@<TRIPOS>BOND
     1     1     2 1
     2     2     3 1
     3     2     5 1
     4     2     6 1
     5     3     4 1
     6     3     7 1
     7     3     8 1
@<TRIPOS>SUBSTRUCTURE
     1 LIG         1 TEMP              0 ****  ****    0 ROOT
"""

FRCMOD = """Remark line goes here
MASS

BOND

ANGLE

DIHE

IMPROPER

NONBON

"""

EXPECTED = {"CL1": ("CL", 17), "C1": ("C", 6), "C2": ("C", 6), "BR1": ("BR", 35),
            "H1": ("H", 1), "H2": ("H", 1), "H3": ("H", 1), "H4": ("H", 1)}


class TestReadAtomTypeElements(unittest.TestCase):

    def test_reads_add_atom_types(self):
        with tempfile.TemporaryDirectory() as tmp:
            leaprc = Path(tmp, "leaprc.test")
            leaprc.write_text('addAtomTypes {\n\t{ "c3"  "C" "sp3" }\n\t{ "cl"  "Cl" "sp3" }\n'
                              '\t{ "br"  "Br" "sp3" }\n\t{ "EP"  "" "sp3" }\n}\n')
            self.assertEqual(_read_atom_type_elements(leaprc), {"c3": "C", "cl": "Cl", "br": "Br"})


@unittest.skipUnless(shutil.which("tleap"), "tleap not available")
class TestStageLeapElements(unittest.TestCase):

    def test_halogen_elements(self):
        with tempfile.TemporaryDirectory() as tmp:
            cwd = Path(tmp)
            mol2 = cwd / "lig.mol2"
            mol2.write_text(MOL2)
            frcmod = cwd / "lig.frcmod"
            frcmod.write_text(FRCMOD)
            StageLeap("Leap", main_input=mol2, cwd=cwd, in_frcmod=frcmod, out_lib=cwd / "lig.lib").execute()

            pdb_elements = {line[12:16].strip(): line[76:78].strip()
                            for line in (cwd / "lig.pdb").read_text().splitlines() if line.startswith("ATOM")}
            self.assertEqual(pdb_elements, {name: element for name, (element, _) in EXPECTED.items()})

            lib_lines = (cwd / "lig.lib").read_text().splitlines()
            start = next(i for i, line in enumerate(lib_lines) if ".unit.atoms table" in line)
            lib_elements = {}
            for line in lib_lines[start + 1:]:
                if line.startswith("!"):
                    break
                fields = line.split()
                lib_elements[fields[0].strip('"')] = int(fields[6])
            self.assertEqual(lib_elements, {name: elmnt for name, (_, elmnt) in EXPECTED.items()})


if __name__ == '__main__':
    unittest.main()
