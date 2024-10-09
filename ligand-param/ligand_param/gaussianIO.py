import os

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
            return False
        
        if os.path.exists(self.filename.strip('.com')+'.log'):
            print(f"File {self.filename.strip('.com')+'.log'} already exists. Exiting.")
            return False

        with open(self.filename, 'w') as f:
            for link in self.links:
                for line in link.generate_block():
                    f.write(f"{line}\n")

        return False

    def print(self):
        for linkno, link in enumerate(self.links):
            if linkno == 1: print("--Link1--")
            link.print()
    
    def add_block(self, block):
        if isinstance(block, GaussianInput):
            self.nlinks += 1
            self.links.append(block)
    
    def get_run_command(self, extension='.com'):
        if extension not in self.filename:
            raise ValueError("Extension does not match filename.")
        return f"g16 < self.filename > {self.filename.strip(extension)}.log"

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




class GaussianReader:
    """ Class to read Gaussian log files """
    def __init__(self, filename):
        """
        Initialize the class with a filename
        """
        self.filename = filename
        return
    
    def read_log(self):
        """ Read the Gaussian log file, and extract information from it.

        This is adapted from ReadGauOutput in Tim's parmutils package.

        Parameters
        ----------
        None

        Returns
        -------
        atn : list
            A list of atomic symbols
        coords : list
            A list of atomic coordinates
        charge : int
            The charge of the molecule
        multiplicity : int
            The multiplicity of the molecule

        """
        atn, coords = [], []
        charge, multiplicity = 0, 1
        readflag = False
        with open(self.filename,'r') as f:
            arc = ""
            for line in f:
                # Document whatever this is
                if readflag:
                    arc += line.strip()
                if '1\\1\\' in line:
                    arc=line.strip()
                    readflag=True
                if '\\\\@' in arc:
                    readflag=False
                    break
            secs = arc.split("\\\\")
        try:
            data = [ sub.split(",") for sub in secs[3].split("\\") ]
            charge = int( data[0][0] )
            multiplicity = int( data[0][1] )
            for i in range( len(data) - 1):
                atn.append(data[i+1][0])
                if len(data[i+1]) == 5:
                    coords.append([float(data[i+1][2]), float(data[i+1][3]), float(data[i+1][4])])
                else:
                    coords.append([float(data[i+1][1]), float(data[i+1][2]), float(data[i+1][3])])
        except:
            raise IOError("Error reading log file")
        
        print(f"Found {len(atn)} atoms.")

        return atn, coords, charge, multiplicity
    



if __name__ == "__main__":
    header = ["%nproc=4", "%mem=2GB"]
    test_coords = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    test_elements = ['C', 'H']
    test = GaussianInput(initial_coordinates=test_coords, elements=test_elements, header=header)


    new_write = GaussianWriter("test.com")
    new_write.add_block(test)
    new_write.add_block(test)
    new_write.write(dry_run=True)

    read_test = GaussianReader("F3KRP.log")
    read_test.read_log()


