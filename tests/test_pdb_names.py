"""Tests for PDB atom-name parsing and pool generation."""

import pytest

from rdkit import Chem

from ligandparam.stages.pdb_names import PDB_Name_Fixer
from ligandparam.stages.smilestopdb import StageSmilesToPDB


def _atom_with_name(name, symbol="C"):
    mol = Chem.MolFromSmiles({"C": "C", "N": "N", "O": "O", "S": "S"}[symbol])
    atom = mol.GetAtomWithIdx(0)
    info = Chem.AtomPDBResidueInfo()
    info.SetName(name)
    atom.SetMonomerInfo(info)
    return atom


PARSERS = pytest.mark.parametrize(
    "parse",
    [PDB_Name_Fixer.get_element_name_and_number, StageSmilesToPDB.get_element_name_and_number],
    ids=["pdb_name_fixer", "smiles_to_pdb"],
)


@PARSERS
@pytest.mark.parametrize("name,symbol", [
    ("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("OXT", "O"), ("SG", "S"),
])
def test_digit_free_names_do_not_raise(parse, name, symbol):
    """Standard PDB atom names carry no digits; int('') used to raise ValueError,
    which aborted any run given a real reference PDB."""
    element_number, parsed_name, number, element = parse(_atom_with_name(name, symbol))
    assert parsed_name == name
    assert number == 0
    assert element == name


@PARSERS
@pytest.mark.parametrize("name,expected_element,expected_number", [
    ("C1", "C", 1),
    ("C12", "C", 12),
    ("CA1", "CA", 1),
])
def test_numbered_names_are_split(parse, name, expected_element, expected_number):
    _, _, number, element = parse(_atom_with_name(name))
    assert (element, number) == (expected_element, expected_number)


@PARSERS
def test_empty_name_falls_back_to_the_element_symbol(parse):
    element_number, name, number, element = parse(_atom_with_name(""))
    assert name == "C0"
    assert number == 0
    assert element == "C"


class TestPadAtomName:
    def test_short_names_are_padded(self):
        assert len(PDB_Name_Fixer.pad_atom_name("C1")) == 4

    def test_four_character_names_are_unchanged(self):
        assert PDB_Name_Fixer.pad_atom_name("C123") == "C123"
