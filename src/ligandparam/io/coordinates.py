import warnings
from functools import lru_cache
from typing import Optional,  Union
import shutil
from pathlib import Path

import numpy as np

import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element, guess_masses


@lru_cache(maxsize=1)
def _ELEMENT_SYMBOLS() -> frozenset:
    """The set of real element symbols, built once from RDKit's periodic table.

    Looked up by symbol rather than probed with ``GetAtomicNumber``, which raises a
    C++ post-condition violation (and prints a long stack trace) for unknown input.
    """
    from rdkit import Chem
    table = Chem.GetPeriodicTable()
    return frozenset(table.GetElementSymbol(z) for z in range(1, 119))


def normalize_element(symbol) -> str:
    """Return ``symbol`` as a canonical element symbol, or "" if it is not one.

    Parameters
    ----------
    symbol : str
        A candidate element symbol, e.g. "N", "cl", "CL".

    Returns
    -------
    str
        The canonical symbol ("N", "Cl"), or an empty string if `symbol` does not
        name a real element.

    Notes
    -----
    MDAnalysis's name-based guesser returns whatever letters it finds rather than
    failing, so an atom named "??" guesses to "?". Validating here keeps such values
    out of the Gaussian input.
    """
    text = str(symbol).strip()
    if not text or not text.isalpha():
        return ""
    candidate = text[0].upper() + text[1:].lower()
    return candidate if candidate in _ELEMENT_SYMBOLS() else ""


def repair_zero_masses(universe: mda.Universe, atol: float = 0.1) -> None:
    """Ensure every atom in ``universe`` has a non-zero mass, in place.

    Mol2 files written by antechamber often carry unusable atom names, so MDAnalysis
    cannot guess masses and leaves them at zero. A zero total mass makes
    ``center_of_mass()`` return NaN, which breaks centering and rotation. The exact
    values do not matter for those operations, so unresolved masses are set to 1.0.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        The universe to repair. Modified in place.
    atol : float, optional
        Absolute tolerance for treating a mass as zero (default 0.1).

    Notes
    -----
    ``AtomGroup.masses`` returns a *copy*, so ``u.atoms.masses[mask] = 1.0`` writes to
    a temporary and silently does nothing. The read-modify-write below is required.
    """
    if np.any(np.isclose(universe.atoms.masses, 0, atol=atol)):
        try:
            universe.guess_TopologyAttrs(to_guess=['elements'], force_guess=['masses'])
        except Exception:
            # Guessing needs names or types to work from and raises when it has
            # neither. That is the case this function exists to survive, so fall
            # through to the 1.0 default below.
            pass

    masses = universe.atoms.masses
    zero = np.isclose(masses, 0, atol=atol)
    if np.any(zero):
        masses[zero] = 1.0
        universe.atoms.masses = masses


class Coordinates:

    def __init__(self, filename: Union[Path, str], filetype: str = 'pdb'):
        """   A class to handle the coordinates of a structure. 
        
        This class is a wrapper around the MDAnalysis Universe class, and provides a simple interface to 
        manipulate the coordinates of a structure.
            
        Parameters
        ----------
        filename : Union[Path, str]
            The filename of the structure to read in
        filetype : str, optional
            The filetype of the structure to read in
        """
        self.filename = Path(filename)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.u = mda.Universe(filename)
        self.original_coords = self.get_coordinates()

        # If the mol2 comes from antechaamber, then the atom names are weird and both rdkit and mda will have trouble
        repair_zero_masses(self.u)

        return

    def get_coordinates(self):
        """ Grabs the coordinates
        
        Parameters
        ----------
        None

        Returns
        -------
        coords : np.array
            The coordinates of the atoms in the structure
        """
        return self.u.atoms.positions

    def get_elements(self):
        """ Grabs the elements

        Every atom is resolved individually: the topology's `elements` attribute is
        used where it holds a value, and the atom name is used to guess the rest.

        Parameters
        ----------
        None

        Returns
        -------
        elements : list of str
            The element symbol of each atom in the structure.

        Raises
        ------
        ValueError
            If any atom's element cannot be determined.

        Notes
        -----
        Reading the `elements` attribute wholesale inside a try/except is
        all-or-nothing: it catches the attribute being absent, but not the far more
        common case of it being *partially* populated. MDAnalysis leaves an empty
        string for any atom it could not resolve -- notably for mol2 files written by
        antechamber, whose GAFF atom types ("na", "cc", "hn") are not element symbols.
        Those empty strings used to flow straight into the Gaussian input, producing
        coordinate lines with no element column at all.
        """
        natoms = len(self.u.atoms)
        raw = getattr(self.u.atoms, "elements", None)
        if raw is None:
            elements = [""] * natoms
        else:
            elements = [normalize_element(e) for e in raw]

        missing = [i for i, e in enumerate(elements) if not e]
        if missing:
            guessed = self._get_elements_from_topology()
            for i in missing:
                if i < len(guessed):
                    elements[i] = normalize_element(guessed[i])

        unresolved = [i for i, e in enumerate(elements) if not e]
        if unresolved:
            names = [str(self.u.atoms[i].name) for i in unresolved]
            raise ValueError(
                f"Could not determine the element for {len(unresolved)} atom(s) in "
                f"{self.filename}: {names}. Writing them out would produce a Gaussian "
                f"input with a blank element column.")

        return elements

    def _get_elements_from_topology(self):
        """ Guesses the elements from the atom names in the topology

        Parameters
        ----------
        None

        Returns
        -------
        elements : list of str
            The guessed element symbol of each atom. An entry may be empty if the
            name could not be interpreted.
        """
        names = self.u.atoms.names
        # MDAnalysis.topology.guessers is deprecated and scheduled for removal in
        # 3.0.0, so prefer the Guesser API and fall back only if it is unavailable.
        try:
            from MDAnalysis.guesser import DefaultGuesser
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                return list(DefaultGuesser(None).guess_atom_element(name) for name in names)
        except ImportError:
            pass

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from MDAnalysis.topology.guessers import guess_types
            return list(guess_types(names))

    def update_coordinates(self, coords, original=False):
        """ Updates the coordinates

        Parameters
        ----------
        coords : np.array
            The new coordinates to update the structure with

        Returns
        -------
        None
        """

        assert np.shape(coords) == np.shape(self.get_coordinates()), "Coordinate dimensions do not match"
        self.u.atoms.positions = coords
        if original:
            self.original_coords = coords
        return

    def rotate(self, alpha=0.0, beta=0.0, gamma=0.0):
        """ Rotate the coordinates around specific axes around the center of mass.

        The rotation is done in the order alpha, beta, gamma, and the rotation is done around the center of mass.
        
        Parameters
        ----------
        alpha : float
            The angle to rotate the structure in the alpha direction (degrees)
        beta : float
            The angle to rotate the structure in the beta direction (degrees)
        use_original : bool, optional
            If True, the rotation will be applied to the new coordinates
        """
        import warnings
        warnings.filterwarnings(
            "ignore")  # There is a deprecation warning that will eventually break this code, but this is something that is broken in MDAnalysis
        import MDAnalysis.transformations

        x, y, z = [1, 0, 0], [0, 1, 0], [0, 0, 1]
        ts = self.u.trajectory.ts

        self.u.atoms.positions = self.original_coords
        com = self.u.atoms.center_of_mass()

        # Apply rotation around the x axis
        rotated = mda.transformations.rotate.rotateby(angle=alpha, direction=x, point=com)(ts)
        self.u.atoms.positions = rotated

        # Apply rotation around the y axis
        rotated = mda.transformations.rotate.rotateby(angle=beta, direction=y, point=com)(ts)
        self.u.atoms.positions = rotated

        # Apply rotation around the z axis
        rotated = mda.transformations.rotate.rotateby(angle=gamma, direction=z, point=com)(ts)
        self.u.atoms.positions = rotated

        return self.get_coordinates()


def SimpleXYZ(file_obj, coordinates):
    """ Write a simple XYZ file with the coordinates. 
    
    Parameters
    ----------
    file_obj : file object
        The file object to write to
    coordinates : np.array
        The coordinates to write to the file
    """
    file_obj.write(f"{len(coordinates)}\n")
    file_obj.write("Generated by ligand_param\n")
    for i, coord in enumerate(coordinates):
        file_obj.write(f"{i + 1} {coord[0]} {coord[1]} {coord[2]}\n")
    return


class Mol2Writer:
    def __init__(self, u, filename=None, selection="all"):
        """ A class to write a mol2 file.
        
        Parameters
        ----------
        u : MDAnalysis Universe
            The universe to write to a mol2 file
        filename : str
            The filename to write to
        """
        self.u = u
        self.filename = Path(filename)
        self.selection = selection
        return

    def _write(self):
        """ Uses MDAnalysis to write the mol2 file. """
        ag = self.u.select_atoms(self.selection)
        ag.write(self.filename)

    def _remove_blank_lines(self):
        """ Remove blank lines from a file.
        
        Parameters
        ----------
        file_path : str
            The path to the file to remove blank lines from
        
        Returns
        -------
        None
        
        """
        if Path(self.filename).exists():
            # Read the file and filter out blank lines
            with open(self.filename, 'r') as file:
                lines = file.readlines()
                non_blank_lines = [line for line in lines if line.strip()]

            # Write the non-blank lines back to the file
            with open(self.filename, 'w') as file:
                file.writelines(non_blank_lines)
        else:
            raise FileNotFoundError(f"File {self.filename} not found.")

    def write(self):
        """ Write the mol2 file. 
        
        This uses the _write method to write the mol2 file, and then removes any blank lines from the file.
        
        Parameters
        ----------
        None
        
        """
        self._write()
        self._remove_blank_lines()
        return


def Remove_PDB_CONECT(filename: Union[Path, str], backup: bool = False):
    """ Removes CONECT lines from a PDB file.

    This script (1) copies the pdb file to a new file (with input added to the filename)
    and (2) removes the CONECT records from the original file.
    
    Parameters
    ----------
    filename : str
        The name of the file to check
        
    Returns
    -------
    None

    """
    fn = Path(filename)
    if backup:
        shutil.copyfile(fn, fn.parent / f"input_{fn.name}")
    with open(filename, 'r') as file:
        lines = file.readlines()
        new_lines = []
        for line in lines:
            if line.strip().startswith("CONECT"):
                continue
            new_lines.append(line)
    with open(filename, 'w') as file:
        file.writelines(new_lines)
    return
