"""Element resolution, and the guard against writing a Gaussian input without elements.

A user reported .com files whose coordinate lines had no element column for most
atoms and a correct symbol for a few. The cause was that the elements topology
attribute can be *partially* populated -- MDAnalysis leaves an empty string for any
atom it cannot resolve, and mol2 files from antechamber carry GAFF atom types
("na", "cc", "hn") rather than element symbols.
"""

import warnings

import numpy as np
import pytest

import MDAnalysis as mda

from ligandparam.io.coordinates import Coordinates
from ligandparam.io.gaussianIO import GaussianInput


MOL2 = """@<TRIPOS>MOLECULE
G0
 4 2 1 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 N9          -1.6180     1.4790    -0.7330 na         1 G0       -0.315900
      2 C8          -1.9110     1.8960     0.5390 cc         1 G0        0.390400
      3 H9          -2.0360     1.7910    -1.5970 hn         1 G0        0.316700
      4 H8          -2.6630     2.6500     0.7250 h5         1 G0        0.068100
@<TRIPOS>BOND
     1    1    2 1
     2    1    3 1
"""


# The reported failure. Uppercase Amber (parm94-style) atom types: MDAnalysis builds
# an elements attribute from the type column here, and only the atoms typed exactly
# "H" resolve -- every other atom comes back as an empty string, including the
# hydrogen typed "H5". The .com that resulted had no element column on 12 of 16 lines.
AMBER_TYPE_MOL2 = """@<TRIPOS>MOLECULE
G0
16 17 1 0 0
SMALL
No Charge or Current Charge
@<TRIPOS>ATOM
   1   N9       -2.5115    2.5603    0.7560   NA 1 G0 -0.3159
   2   C8       -2.8035    2.9932    2.0700   CR 1 G0  0.3904
   3   N7       -2.0525    2.3842    3.0070   NB 1 G0 -0.5830
   4   C5       -1.2285    1.5113    2.2940   CB 1 G0  0.1758
   5   C6       -0.2035    0.5963    2.7520    C 1 G0  0.6520
   6   O6        0.1725    0.3743    3.9110    O 1 G0 -0.5965
   7   N1        0.4245   -0.1347    1.7050   NA 1 G0 -0.5094
   8   C2        0.1035    0.0482    0.3450   CA 1 G0  0.6838
   9   N2        0.8905   -0.6428   -0.6030   N2 1 G0 -0.8928
  10   N3       -0.8675    0.8782   -0.1090   NC 1 G0 -0.6751
  11   C4       -1.5005    1.6083    0.8810   CB 1 G0  0.1088
  12   H9       -2.9275    2.8653   -0.0830    H 1 G0  0.3167
  13   H8       -3.5715    3.7493    2.2630   H5 1 G0  0.0681
  14  H21        1.1825   -1.5568   -0.3360    H 1 G0  0.4213
  15  H22        0.5645   -0.5768   -1.5440    H 1 G0  0.4213
  16   H1        1.1805   -0.7318    1.9630    H 1 G0  0.3345

@<TRIPOS>BOND
    1     1     2  1
    2     1    11  1
    3     1    12  1
    4     2     3  2
    5     2    13  1
    6     3     4  1
    7     4     5  1
    8     4    11  2
    9     5     6  2
   10     5     7  1
   11     7     8  1
   12     7    16  1
   13     8     9  1
   14     8    10  2
   15     9    14  1
   16     9    15  1
   17    10    11  1
"""

EXPECTED_G0_ELEMENTS = ["N", "C", "N", "C", "C", "O", "N", "C", "N", "N", "C",
                        "H", "H", "H", "H", "H"]


@pytest.fixture
def antechamber_mol2(tmp_path):
    path = tmp_path / "lig.mol2"
    path.write_text(MOL2)
    return path


@pytest.fixture
def amber_type_mol2(tmp_path):
    path = tmp_path / "G0.mol2"
    path.write_text(AMBER_TYPE_MOL2)
    return path


def test_uppercase_amber_types_resolve(amber_type_mol2):
    """Regression: uppercase Amber types left 12 of 16 elements blank in the .com."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        elements = Coordinates(amber_type_mol2, filetype="mol2").get_elements()
    assert list(elements) == EXPECTED_G0_ELEMENTS


def test_amber_types_are_not_read_as_metals(amber_type_mol2):
    """NA/CR/NB/CA are valid symbols for sodium/chromium/niobium/calcium; here they
    are nitrogen and carbon, and the atom name is the reliable source."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        elements = Coordinates(amber_type_mol2, filetype="mol2").get_elements()
    assert not ({"Na", "Cr", "Nb", "Ca"} & set(elements))


def test_gaussian_block_from_amber_types_has_every_element(amber_type_mol2):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        coords = Coordinates(amber_type_mol2, filetype="mol2")
        block = GaussianInput(
            command="#P HF/6-31G* OPT(CalcFC)",
            initial_coordinates=coords.get_coordinates(),
            elements=coords.get_elements(), charge=0, header=["%NPROC=8"])
    atom_lines = [ln for ln in block.generate_block()
                  if ln.startswith("     ") and ln.strip()]
    assert len(atom_lines) == 16
    for line, expected in zip(atom_lines, EXPECTED_G0_ELEMENTS):
        assert line.split()[0] == expected


def test_gaff_types_are_not_mistaken_for_elements(antechamber_mol2):
    """GAFF types like 'na' must not resolve to sodium, nor to a blank."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        elements = Coordinates(antechamber_mol2, filetype="mol2").get_elements()
    assert list(elements) == ["N", "C", "H", "H"]


def test_no_element_is_blank(antechamber_mol2):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        elements = Coordinates(antechamber_mol2, filetype="mol2").get_elements()
    assert all(str(e).strip() for e in elements)


def test_partially_populated_elements_are_filled_in(tmp_path, monkeypatch):
    """The failure mode from the report: some atoms resolved, others left empty.

    Reading the attribute wholesale inside a try/except cannot catch this, because
    no exception is raised.
    """
    path = tmp_path / "lig.mol2"
    path.write_text(MOL2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        coords = Coordinates(path, filetype="mol2")
    # Simulate a topology that resolved only the 'hn' hydrogen.
    coords.u.add_TopologyAttr("elements", ["", "", "H", ""])
    assert list(coords.get_elements()) == ["N", "C", "H", "H"]


def test_unresolvable_elements_raise_with_the_atom_names(tmp_path):
    path = tmp_path / "lig.mol2"
    path.write_text(MOL2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        coords = Coordinates(path, filetype="mol2")
    coords.u.add_TopologyAttr("elements", ["", "", "", ""])
    coords.u.add_TopologyAttr("names", ["??", "??", "??", "??"])
    with pytest.raises(ValueError, match="Could not determine the element"):
        coords.get_elements()


class TestGaussianInputGuard:
    COORDS = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])

    def test_blank_element_is_refused(self):
        block = GaussianInput(command="#P HF/6-31G*", initial_coordinates=self.COORDS,
                              elements=["C", ""], charge=0, header=["%NPROC=8"])
        with pytest.raises(ValueError, match="Missing element symbol"):
            block.generate_block()

    def test_whitespace_element_is_refused(self):
        block = GaussianInput(command="#P HF/6-31G*", initial_coordinates=self.COORDS,
                              elements=["C", "  "], charge=0, header=["%NPROC=8"])
        with pytest.raises(ValueError, match="Missing element symbol"):
            block.generate_block()

    def test_complete_elements_are_written(self):
        block = GaussianInput(command="#P HF/6-31G*", initial_coordinates=self.COORDS,
                              elements=["C", "H"], charge=0, header=["%NPROC=8"])
        atom_lines = [ln for ln in block.generate_block() if ln.strip().startswith(("C ", "H "))]
        assert len(atom_lines) == 2

    def test_geometry_free_block_is_still_allowed(self):
        """The Link1 RESP block carries no geometry and must keep working."""
        block = GaussianInput(command="#P HF/6-31G* GEOM(AllCheck) Guess(Read)",
                              charge=0, header=["%NPROC=8"])
        assert block.generate_block()
