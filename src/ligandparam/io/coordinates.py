import warnings
from typing import Optional,  Union
import shutil
from pathlib import Path

import numpy as np

import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element, guess_masses


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
        if np.any(np.isclose(self.u.atoms.masses, 0, atol=0.1)):
            self.u.guess_TopologyAttrs(to_guess=['elements'], force_guess=['masses'])
        # We tried to get correct masses but may have failed in the process. Lack of masses will fail
        # MDAnalysis's center_of_mass(), so just set them to 1.0, since the exact values are not important
        self.u.atoms.masses[np.isclose(self.u.atoms.masses, 0, atol=0.1)] = 1.0

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

        Parameters
        ----------
        None

        Returns
        -------
        elements : list
            The elements of the atoms in the structure
        """
        try:
            return [atom.element for atom in self.u.atoms]
        except:
            return self._get_elements_from_topology()

    def _get_elements_from_topology(self):
        """ Grabs the elements from the topology
        
        Parameters
        ----------
        None

        Returns
        -------
        elements : list
            The elements of the atoms in the structure
        """
        from MDAnalysis.topology.guessers import guess_types
        elements = guess_types(self.u.atoms.names)
        return elements

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
