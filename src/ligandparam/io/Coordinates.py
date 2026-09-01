import warnings
from typing import Optional,  Union
import shutil
from pathlib import Path

import numpy as np

import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element, guess_masses


class Coordinates:
    """Thin MDAnalysis wrapper for reading and transforming structure coordinates."""

    def __init__(self, filename: Union[Path, str], filetype: str = 'pdb'):
        """
        Load a structure and sanitize masses for center-of-mass operations.

        Parameters
        ----------
        filename : Union[Path, str]
            Path to the structure file to read.
        filetype : str, optional
            File type hint (default: ``'pdb'``).
        """
        self.filename = Path(filename)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.u = mda.Universe(filename)
        self.original_coords = np.array(self.get_coordinates(), dtype=float, copy=True)

        # If the mol2 comes from antechaamber, then the atom names are weird and both rdkit and mda will have trouble
        if np.any(np.isclose(self.u.atoms.masses, 0, atol=0.1)):
            self.u.guess_TopologyAttrs(to_guess=['elements'], force_guess=['masses'])
        # We tried to get correct masses but may have failed in the process. Lack of masses will fail
        # MDAnalysis's center_of_mass(), so just set them to 1.0, since the exact values are not important
        self.u.atoms.masses[np.isclose(self.u.atoms.masses, 0, atol=0.1)] = 1.0

        return

    def get_coordinates(self):
        """Return the current atomic coordinates.

        Returns
        -------
        np.ndarray
            Coordinates of the atoms in the structure.
        """
        return self.u.atoms.positions

    def get_elements(self):
        """Return atomic element symbols.

        Returns
        -------
        list
            Element symbols for each atom.
        """
        try:
            return [atom.element for atom in self.u.atoms]
        except (AttributeError, TypeError, ValueError) as exc:
            import warnings
            warnings.warn(
                f"Could not read atom.element ({type(exc).__name__}: {exc}); "
                "guessing elements from atom names",
                stacklevel=2,
            )
            return self._get_elements_from_topology()

    def _get_elements_from_topology(self):
        """Guess element symbols from atom names in the topology.

        Returns
        -------
        list
            Guessed element symbols for each atom.
        """
        from MDAnalysis.topology.guessers import guess_types
        elements = guess_types(self.u.atoms.names)
        return elements

    def update_coordinates(self, coords, original=False):
        """Replace the current atomic coordinates.

        Parameters
        ----------
        coords : np.ndarray
            New coordinates with the same shape as the current positions.
        original : bool, optional
            If True, also update the stored original coordinates used by
            :meth:`rotate`.
        """
        assert np.shape(coords) == np.shape(self.get_coordinates()), "Coordinate dimensions do not match"
        self.u.atoms.positions = coords
        if original:
            self.original_coords = coords
        return

    def rotate(self, alpha=0.0, beta=0.0, gamma=0.0):
        """Rotate coordinates about the center of mass using Euler angles.

        Rotations are applied in order alpha (x), beta (y), gamma (z), matching
        the previous MDAnalysis ``rotateby`` sequence. Angles are in degrees.

        Parameters
        ----------
        alpha : float
            Rotation about the x-axis (degrees).
        beta : float
            Rotation about the y-axis (degrees).
        gamma : float
            Rotation about the z-axis (degrees).

        Returns
        -------
        np.ndarray
            Rotated coordinates with shape ``(n_atoms, 3)``.
        """
        coords = np.asarray(self.original_coords, dtype=float)
        # COM from original geometry (masses already sanitized in __init__)
        masses = self.u.atoms.masses
        com = np.average(coords, axis=0, weights=masses)

        a, b, g = np.deg2rad([alpha, beta, gamma])
        ca, sa = np.cos(a), np.sin(a)
        cb, sb = np.cos(b), np.sin(b)
        cg, sg = np.cos(g), np.sin(g)

        # Intrinsic/extrinsic composition matching sequential Rx, Ry, Rz on positions
        rx = np.array([[1.0, 0.0, 0.0],
                       [0.0, ca, -sa],
                       [0.0, sa, ca]])
        ry = np.array([[cb, 0.0, sb],
                       [0.0, 1.0, 0.0],
                       [-sb, 0.0, cb]])
        rz = np.array([[cg, -sg, 0.0],
                       [sg, cg, 0.0],
                       [0.0, 0.0, 1.0]])
        rotation = rz @ ry @ rx

        rotated = (coords - com) @ rotation.T + com
        self.u.atoms.positions = rotated
        return rotated

    def rotate_matrix(self, rotation: np.ndarray) -> np.ndarray:
        """Rotate the original coordinates using an explicit rotation matrix.

        The rotation is applied about the mass-weighted center of mass. This
        is the path used by quaternion SO(3) orientation protocols.

        Parameters
        ----------
        rotation : np.ndarray, shape (3, 3)
            Proper orthogonal rotation matrix.

        Returns
        -------
        np.ndarray
            Rotated coordinates with shape ``(n_atoms, 3)``.

        Raises
        ------
        ValueError
            If ``rotation`` is not a valid proper rotation matrix.
        """
        rotation = np.asarray(rotation, dtype=float)
        if rotation.shape != (3, 3):
            raise ValueError(f"Expected rotation shape (3, 3), got {rotation.shape}")
        if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-8):
            raise ValueError("Rotation matrix must be orthogonal")
        if not np.isclose(np.linalg.det(rotation), 1.0, atol=1e-8):
            raise ValueError("Rotation matrix must have determinant +1")

        coords = np.asarray(self.original_coords, dtype=float)
        com = np.average(coords, axis=0, weights=self.u.atoms.masses)
        rotated = (coords - com) @ rotation.T + com
        self.u.atoms.positions = rotated
        return rotated


def SimpleXYZ(file_obj, coordinates):
    """Write coordinates to a simple XYZ trajectory frame.

    Parameters
    ----------
    file_obj : file object
        Open file handle to write to.
    coordinates : np.ndarray
        Coordinates to write.
    """
    file_obj.write(f"{len(coordinates)}\n")
    file_obj.write("Generated by ligand_param\n")
    for i, coord in enumerate(coordinates):
        file_obj.write(f"{i + 1} {coord[0]} {coord[1]} {coord[2]}\n")
    return


class Mol2Writer:
    """Write an MDAnalysis Universe selection to a mol2 file."""

    def __init__(self, u, filename=None, selection="all"):
        """
        Parameters
        ----------
        u : MDAnalysis.Universe
            Universe to write.
        filename : str, optional
            Output mol2 path.
        selection : str, optional
            Atom selection string (default: ``'all'``).
        """
        self.u = u
        self.filename = Path(filename)
        self.selection = selection
        return

    def _write(self):
        """Write the selected atoms to mol2 via MDAnalysis."""
        ag = self.u.select_atoms(self.selection)
        ag.write(self.filename)

    def _remove_blank_lines(self):
        """Remove blank lines from the written mol2 file.

        Raises
        ------
        FileNotFoundError
            If the output file does not exist.
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
        """Write the mol2 file and strip blank lines that confuse antechamber."""
        self._write()
        self._remove_blank_lines()
        return


def Remove_PDB_CONECT(filename: Union[Path, str], backup: bool = False):
    """Remove CONECT records from a PDB file in place.

    Optionally copies the original to ``input_<filename>`` before editing.

    Parameters
    ----------
    filename : Union[Path, str]
        PDB file to clean.
    backup : bool, optional
        If True, save a copy before removing CONECT lines.
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
