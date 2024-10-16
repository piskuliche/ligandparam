import MDAnalysis as mda
import numpy as np

class Coordinates:
    """ A class to handle the coordinates of a structure. 
    
    This class is a wrapper around the MDAnalysis Universe class, and provides a simple interface to 
    manipulate the coordinates of a structure.
    
    """
    def __init__(self, filename, filetype='pdb'):
        """ Initialize the coordinates object.
        
        Parameters
        ----------
        filename : str
            The filename of the structure to read in
        filetype : str, optional
            The filetype of the structure to read in
        """
        self.filename = filename
        self.u = mda.Universe(filename)
        self.original_coords = self.get_coordinates()
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

        return [atom.element for atom in self.u.atoms]
    
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
        warnings.filterwarnings("ignore") # There is a deprecation warning that will eventually break this code, but this is something that is broken in MDAnalysis
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
        file_obj.write(f"{i+1} {coord[0]} {coord[1]} {coord[2]}\n")
    return



class Mol2Writer:
    """ A class to write a mol2 file. """
    def __init__(self, u, filename=None, selection="all"):
        """ Initialize the Mol2Writer class.
        
        Parameters
        ----------
        u : MDAnalysis Universe
            The universe to write to a mol2 file
        filename : str
            The filename to write to
        """
        self.u = u
        if not Path(filename).exists():
            raise FileNotFoundError(f"File {filename} not found.")
        self.filename = filename
        return
    
    def _write(self):
        """ Uses MDAnalysis to write the mol2 file. """
        ag = self.u.select_atoms(selection)
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

