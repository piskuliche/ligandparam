import MDAnalysis as mda
import numpy as np

class Coordinates:
    def __init__(self, filename, filetype='pdb'):
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
    
    def update_coordinateS(self, coords, original=False):
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

    def rotate(self, alpha=0.0, beta=0.0, use_original=False):
        """ Rotate the coordinates 
        
        Parameters
        ----------
        alpha : float
            The angle to rotate the structure in the alpha direction (degrees)
        beta : float
            The angle to rotate the structure in the beta direction (degrees)
        use_original : bool, optional
            If True, the rotation will be applied to the new coordinates
        """
        output_coords = None
        if use_original is True:
            output_coords = self.original_coords
        return output_coords

def Rotate(coords, alpha=0.0, beta=0.0):
    pass
