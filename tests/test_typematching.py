"""Tests for the mol2 atom-name rewriting used by StageMatchAtomNames."""

import pytest

from ligandparam.stages.typematching import rename_mol2_atom_names


MOL2 = """@<TRIPOS>MOLECULE
LIG
 3 2 1 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 C12       1.0000    2.0000    3.0000 c3        1 LIG      -0.1000
      2 C1        4.0000    5.0000    6.0000 c3        1 LIG       0.0500
      3 H1        7.0000    8.0000    9.0000 hc        1 LIG       0.0500
@<TRIPOS>BOND
     1    1    2 1
     2    1    3 1
"""


def _atom_lines(lines):
    out, seen = [], False
    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            seen = True
            continue
        if seen:
            if line.startswith("@<TRIPOS>"):
                break
            out.append(line)
    return out


def _names(lines):
    return [line.split()[1] for line in _atom_lines(lines)]


def test_shorter_replacement_does_not_leave_a_fragment():
    """"C12" -> "CA" used to yield "CA2": the splice was sized by the new name."""
    out, exhausted = rename_mol2_atom_names(MOL2.splitlines(keepends=True), ["CA", "CB", "HA"])
    assert not exhausted
    assert _names(out) == ["CA", "CB", "HA"]


def test_longer_replacement_does_not_eat_the_coordinate_column():
    out, exhausted = rename_mol2_atom_names(MOL2.splitlines(keepends=True), ["CARBON1", "CB", "HA"])
    assert not exhausted
    assert _names(out) == ["CARBON1", "CB", "HA"]
    # The x coordinate of the first atom must survive intact.
    assert _atom_lines(out)[0].split()[2] == "1.0000"


def test_all_coordinates_are_preserved():
    original = MOL2.splitlines(keepends=True)
    out, _ = rename_mol2_atom_names(original, ["CA", "CB", "HA"])
    for before, after in zip(_atom_lines(original), _atom_lines(out)):
        assert before.split()[2:] == after.split()[2:]


def test_bond_section_is_untouched():
    """The regex matches bond records too; renaming must stop at the section header."""
    original = MOL2.splitlines(keepends=True)
    out, exhausted = rename_mol2_atom_names(original, ["CA", "CB", "HA"])
    assert not exhausted, "bond records were consuming names"
    assert out[out.index("@<TRIPOS>BOND\n"):] == original[original.index("@<TRIPOS>BOND\n"):]


def test_header_is_preserved():
    original = MOL2.splitlines(keepends=True)
    out, _ = rename_mol2_atom_names(original, ["CA", "CB", "HA"])
    stop = original.index("@<TRIPOS>ATOM\n") + 1
    assert out[:stop] == original[:stop]


def test_too_few_names_is_reported():
    out, exhausted = rename_mol2_atom_names(MOL2.splitlines(keepends=True), ["CA"])
    assert exhausted
    # The atoms that could be renamed still were; the rest keep their original names.
    assert _names(out) == ["CA", "C1", "H1"]


def test_extra_names_are_ignored():
    out, exhausted = rename_mol2_atom_names(
        MOL2.splitlines(keepends=True), ["CA", "CB", "HA", "HB", "HC"])
    assert not exhausted
    assert _names(out) == ["CA", "CB", "HA"]
