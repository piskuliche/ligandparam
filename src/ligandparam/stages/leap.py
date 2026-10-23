
"""
StageLeap module
----------------
This module provides the StageLeap class for generating leap library files for molecular systems.
"""

from typing import Optional,  Union, Any

from pathlib import Path
import os
import re
import shutil
import warnings

import MDAnalysis as mda

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.io.leapIO import LeapWriter
from ligandparam.interfaces import Leap
from ligandparam.utils import find_word_and_get_line

# GAFF/GAFF2 types that loadmol2 assigns the wrong element (it uses the first letter of the type: cl -> C, br -> B).
# Used when a leaprc's addAtomTypes table can't be read.
_FALLBACK_TYPE_ELEMENTS = {"cl": "Cl", "br": "Br"}
# An addAtomTypes entry, e.g. { "cl"  "Cl" "sp3" }
_ADD_ATOM_TYPE_RE = re.compile(r'\{\s*"([^"]+)"\s+"([^"]*)"\s+"[^"]*"\s*\}')


def _find_leaprc(leaprc: Union[Path, str]) -> Optional[Path]:
    """
    Find a leaprc file the way tleap would, either as a path or in the AmberTools leap/cmd directory.

    Parameters
    ----------
    leaprc : Union[Path, str]
        A path to a leaprc file, or its name (e.g. 'leaprc.gaff2').

    Returns
    -------
    Optional[Path]
        The path to the leaprc file, or None if it can't be found.
    """
    if Path(leaprc).is_file():
        return Path(leaprc)
    cmd_dirs = []
    if amberhome := os.environ.get("AMBERHOME"):
        cmd_dirs.append(Path(amberhome, "dat", "leap", "cmd"))
    if tleap := shutil.which("tleap"):
        cmd_dirs.append(Path(tleap).resolve().parents[1] / "dat" / "leap" / "cmd")
    for cmd_dir in cmd_dirs:
        if Path(cmd_dir, leaprc).is_file():
            return Path(cmd_dir, leaprc)
    return None


def _read_atom_type_elements(leaprc: Path) -> dict[str, str]:
    """
    Read the atom type to element mapping from the addAtomTypes entries of a leaprc file.

    Parameters
    ----------
    leaprc : Path
        Path to the leaprc file.

    Returns
    -------
    dict[str, str]
        Mapping of atom type to element symbol. Types without an element (extra points) are skipped.
    """
    return {atom_type: element for atom_type, element in _ADD_ATOM_TYPE_RE.findall(Path(leaprc).read_text())
            if element}


def _write_pdb_elements(pdb: Path, names: list[str], elements: list[Optional[str]]) -> bool:
    """
    Write the element column (columns 77-78) of the ATOM/HETATM records of a pdb file.

    Parameters
    ----------
    pdb : Path
        Path to the pdb file, which is modified in place.
    names : list[str]
        Expected atom names, in pdb order.
    elements : list[Optional[str]]
        Element symbol for each atom. Atoms with None keep their current element column.

    Returns
    -------
    bool
        True if the pdb was updated, False if its atoms don't match `names`, in which case it is left unchanged.
    """
    lines = Path(pdb).read_text().splitlines(keepends=True)
    atom_lines = [i for i, line in enumerate(lines) if line.startswith(("ATOM  ", "HETATM"))]
    if len(atom_lines) != len(names) or any(lines[i][12:16].strip() != name[:4] for i, name in zip(atom_lines, names)):
        return False
    for i, element in zip(atom_lines, elements):
        if element:
            line = lines[i].rstrip("\n")
            lines[i] = f"{line[:76]:<76}{element.upper():>2}{line[78:]}\n"
    Path(pdb).write_text("".join(lines))
    return True


class StageLeap(AbstractStage):
    """
    Generate leap library files for molecular systems using AmberTools.

    LEaP's loadmol2 guesses elements from the atom types and gets two-letter elements wrong (GAFF ``cl``/``br`` become
    C/B), so the elements in the output lib and pdb are set from the leaprc files' addAtomTypes tables instead.

    Parameters
    ----------
    stage_name : str
        The name of the stage.
    main_input : Union[Path, str]
        Path to the input mol2 file.
    cwd : Union[Path, str]
        Current working directory.
    in_frcmod : str
        Path to the input frcmod file.
    out_lib : str
        Path to the output library file.
    molname : str, optional
        Molecule name (default: 'MOL').
    leaprc : list, optional
        List of leaprc configuration files (default: ['leaprc.gaff2']).

    Attributes
    ----------
    in_mol2 : Path
        Path to the input mol2 file.
    in_frcmod : str
        Path to the input frcmod file.
    out_lib : Path
        Path to the output library file.
    molname : str
        Molecule name.
    leaprc : list
        List of leaprc configuration files.
    out_pdb : Path
        Path to the output pdb file.
    """
    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        """
        Initialize the StageLeap class.

        Parameters
        ----------
        stage_name : str
            The name of the stage.
        main_input : Union[Path, str]
            Path to the input mol2 file.
        cwd : Union[Path, str]
            Current working directory.
        *args
            Additional positional arguments.
        **kwargs
            Additional keyword arguments.
        """
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_mol2 = Path(main_input)
        self.add_required(self.in_mol2)
        self.in_frcmod = kwargs["in_frcmod"]
        self.add_required(self.in_frcmod)
        self.out_lib = Path(kwargs["out_lib"])
        self.molname = kwargs.get("molname", "MOL")
        self.leaprc = kwargs.get("leaprc", ["leaprc.gaff2"])

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """
        Append a stage to the current stage.

        Parameters
        ----------
        stage : AbstractStage
            The stage to append.

        Returns
        -------
        AbstractStage
            The appended stage.
        """
        return stage



    def execute(self, dry_run=False, nproc: Optional[int]=None, mem: Optional[int]=None) -> Any:
        """
        Setup and execute the leap library file generation.

        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would be run.
        nproc : int, optional
            Number of processors to use.
        mem : int, optional
            Amount of memory to use (in GB).

        Returns
        -------
        None
        """
        # Generate the leap input file
        leapgen = LeapWriter("param")
        # Add the leaprc files
        for rc in self.leaprc:
            leapgen.add_leaprc(rc)

        # Try not to overwrite an existing pdb file
        self.out_pdb = Path(self.out_lib.parent, f"{self.out_lib.stem}.pdb")
        if self.out_pdb.is_file():
            self.out_pdb = Path(self.out_lib.parent, f"tleap_{self.out_lib.stem}.pdb")

        # Add the leap commands
        leapgen.add_line(f"loadamberparams {self.in_frcmod.name}")
        leapgen.add_line(f"{self.molname} = loadmol2 {self.in_mol2.name}")
        leapgen.add_line(f'set {self.molname}.1 name "{self.molname}"')
        # Like the residue rename above, this assumes a single-residue mol2.
        names, elements = self._atom_elements()
        for i, element in enumerate(elements, start=1):
            if element:
                leapgen.add_line(f'set {self.molname}.1.{i} element "{element}"')
        leapgen.add_line(f"saveOff {self.molname} {self.out_lib}")
        leapgen.add_line(f"savePDB {self.molname} {self.out_pdb}")
        leapgen.add_line("quit")
        # Write the leap input file
        leapgen.write(self.cwd / "tleap.param.in")
        leap_log = Path(self.cwd, "leap.log")
        leap_log.unlink(missing_ok=True)
        if self.out_lib.is_file():
            shutil.move(self.out_lib, self.out_lib.with_name(f"backup_{self.out_lib.name}"))
        # Call the leap program
        leap = Leap(cwd=self.cwd, logger=self.logger)
        leap.call(f=self.cwd / "tleap.param.in", dry_run=dry_run)

        if lines := find_word_and_get_line(leap_log, "Warning!"):
            self.logger.warning(f"Warning! found in {leap_log}\n{lines}")

        # Not every tleap version writes the pdb element column, and without it readers guess Cl/Br from the atom
        # name as C/B, so write it explicitly.
        out_pdb = self.out_pdb if self.out_pdb.is_absolute() else self.cwd / self.out_pdb
        if not dry_run and any(elements) and not _write_pdb_elements(out_pdb, names, elements):
            self.logger.warning(f"Atoms in {out_pdb} don't match {self.in_mol2.name}; its elements were not set.")

        return

    def _atom_elements(self) -> tuple[list[str], list[Optional[str]]]:
        """
        Get the atom names and elements of the input mol2, using the atom type to element mapping of the leaprc files.

        Returns
        -------
        tuple[list[str], list[Optional[str]]]
            The atom names, and the element of each atom (None if its type has no known element).
        """
        mol2 = self.cwd / self.in_mol2.name
        if not mol2.is_file():
            # Nothing to read yet, e.g. in a dry run.
            return [], []

        type_elements = dict(_FALLBACK_TYPE_ELEMENTS)
        for rc in self.leaprc:
            if (rc_path := _find_leaprc(rc)) is None:
                self.logger.warning(f"Could not find {rc} to read its atom type elements.")
                continue
            type_elements.update(_read_atom_type_elements(rc_path))

        # Supress the inevitable mol2 file warnings.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            atoms = mda.Universe(mol2, format="mol2").atoms
        elements = [type_elements.get(atom_type) for atom_type in atoms.types]
        if unknown := sorted({t for t, el in zip(atoms.types, elements) if el is None}):
            self.logger.warning(f"No element known for atom types {unknown}; LEaP will guess their elements.")
        return list(atoms.names), elements

    def _clean(self):
        """
        Clean the files generated during the stage.
        """
        raise NotImplementedError("clean method not implemented")
