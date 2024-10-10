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
        import MDAnalysis.transformations

        x, y, z = [1, 0, 0], [0, 1, 0], [0, 0, 1]
        ts = self.u.trajectory.ts

        self.u.atoms.positions = self.original_coords
        com = self.u.atoms.center_of_mass()

        rotated = mda.transformations.rotate.rotateby(angle=alpha, direction=x, point=com)(ts)
        self.u.atoms.positions = rotated

        rotated = mda.transformations.rotate.rotateby(angle=beta, direction=y, point=com)(ts)
        self.u.atoms.positions = rotated

        rotated = mda.transformations.rotate.rotateby(angle=gamma, direction=z, point=com)(ts)
        self.u.atoms.positions = rotated

        print(np.shape(self.get_coordinates()))
        return self.get_coordinates()


