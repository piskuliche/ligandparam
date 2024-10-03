import numpy as np


class GaussianWriter:
    """ Class to write Gaussian input files """
    def __init__(self, filename):
        """ Initialize the class with a filename
        
        Parameters
        ----------
        filename : str
            The name of the file to write
        
        Returns
        -------
        None
        
        """
        self.filename = filename
        self.nlinks = 0
        self.links = []
        return
    
    def write(self, dry_run=False):
        """ Write the Gaussian input file 
        
        Parameters
        ----------
        dry_run : bool, optional
            If True, the file will not be written to disk
        
        Returns
        -------
        None
        
        """
        if dry_run: 
            self.print()
            return True

        with open(self.filename, 'w') as f:
            for link in self.links:
                for line in link.generate_block():
                    f.write(line)

        return True

    def print(self):
        for linkno, link in enumerate(self.links):
            if linkno == 1: print("--Link1--")
            link.print()
    
    def add_block(self, block):
        if isinstance(block, GaussianInput):
            self.nlinks += 1
            self.links.append(block)

class GaussianInput:
    """ Class to represent a Gaussian LINK1 block """
    def __init__(self, command="# HF/6-31G* OPT", elements=None, initial_coordinates=None, charge=0, multiplicity=1, header=None):


        if initial_coordinates is not None:
            assert elements is not None, "Elements must be specified if coordinates are provided"
            assert np.shape(initial_coordinates)[0] == len(elements), "Number of elements and coordinates do not match"

        self.command = command
        self.elements = elements
        self.coords = initial_coordinates
        self.charge = charge
        self.multiplicity = multiplicity
        self.header = header

        return
    
    def __str__(self):
        print(self)
        return
    
    def generate_block(self):
        """ Print the Gaussian input block with the information stored in the class
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        lines = []
        if self.header:
            for line in self.header:
                lines.append(line)
        lines.append(f"{self.command}\n")
        lines.append("Gaussian Calculation\n")
        lines.append(f"{self.charge} {self.multiplicity}")
        if self.elements is not None:
            for i, element in enumerate(self.elements):
                lines.append(f"     {element} {self.coords[i][0]: >8.5f} {self.coords[i][1]: >8.5f} {self.coords[i][2]: >8.5f} ")
        lines.append("\n")

        return lines
    
    def print(self):
        for line in self.generate_block():
            print(line)




if __name__ == "__main__":
    header = ["%nproc=4", "%mem=2GB"]
    test_coords = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    test_elements = ['C', 'H']
    test = GaussianInput(initial_coordinates=test_coords, elements=test_elements, header=header)


    new_write = GaussianWriter("test.com")
    new_write.add_block(test)
    new_write.add_block(test)
    new_write.write(dry_run=True)
