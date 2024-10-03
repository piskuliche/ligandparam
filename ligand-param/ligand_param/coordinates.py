import MDAnalysis as mda

class Coordinates:
    def __init__(self, filename, filetype='pdb'):
        self.filename = filename
        self.u = mda.Universe(filename)
        return

    def get_coordinates(self):
        return self.u.atoms.positions

    def get_elements(self):
        return [atom.element for atom in self.u.atoms]
    
    def rotate(self, axis, angle):
        pass
