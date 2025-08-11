import os

import numpy as np
from pathlib import Path


class GaussianWriter:
    def \
            __init__(self, filename):
        """ Class for writing Gaussian input files

        The filename selected will be the name of the Gaussian input file that is written to 
        disk. The class also initializes the number of links to zero, and creates an empty list
        to store the GaussianInput objects.
        
        Parameters
        ----------
        filename : str
            The name of the file to write
        
        Returns
        -------
        None
        
        """
        self.filename = filename
        self.out_dir = Path(filename).parent.mkdir(exist_ok=True)
        self.nlinks = 0
        self.links = []
        return
    
    def write(self, dry_run=False):
        """ Write the Gaussian input file to a file

        This function writes the Gaussian input file to disk. If dry_run is set to True, the file will not be written
        to disk, but will be printed to the screen instead. This is useful for debugging. The file is written in the
        Gaussian input file format, with the LINK1 blocks separated by the --Link1-- delimiter.
        
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
            return 

        with open(self.filename, 'w') as f:
            for i, link in enumerate(self.links):
                if i > 0:
                    f.write("--Link1--\n")
                for line in link.generate_block():
                    f.write(f"{line}\n")

        return 

    def print(self):
        """ Print the Gaussian input file to the screen

        This function prints the Gaussian input file to the screen. This is useful for debugging, or for
        checking the contents of the file before writing it to disk.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        for linkno, link in enumerate(self.links):
            if linkno == 1: print("--Link1--")
            link.print()
    
    def add_block(self, block):
        """ Add a GaussianInput block to the GaussianWriter
        
        This function adds a GaussianInput block to the GaussianWriter. This block will be written to the
        Gaussian input file when the write() function is called.
        
        Parameters
        ----------
        block : GaussianInput
            The GaussianInput block to add to the GaussianWriter
        
        Returns
        -------
        None
        
        """
        if isinstance(block, GaussianInput):
            self.nlinks += 1
            self.links.append(block)
    
    def get_run_command(self, extension='.com'):
        """ Get the command to run the Gaussian input file
        
        This function returns a string that can be used to run the Gaussian input file. This is useful for
        debugging, or for running the Gaussian input file from a script.
        
        Parameters
        ----------
        extension : str, optional
            The extension of the Gaussian input file
        
        Returns
        -------
        str
            The command to run the Gaussian input file
        
        """
        if extension not in self.filename:
            raise ValueError("Extension does not match filename.")
        return f"g16 < self.filename > {self.filename.strip(extension)}.log"

class GaussianInput:
    def __init__(self, command="# HF/6-31G* OPT", elements=None, initial_coordinates=None, charge=0, multiplicity=1, header=None):
        """ Initialize a Gaussian block with the specified parameters.

        Parameters
        ----------
        command : str, optional
            The command for the Gaussian calculation
        elements : list, optional
            A list of atomic symbols
        initial_coordinates : np.array, optional
            A numpy array of atomic coordinates
        charge : int, optional
            The charge of the molecule
        multiplicity : int, optional
            The multiplicity of the molecule
        header : list, optional
            A list of strings to be included in the header of the Gaussian input file

        """

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
        """ Generates the gaussina input block as a list of strings

        This function generates the Gaussian input block as a list of strings. This is useful for writing
        the block to a file, or printing it to the screen.
        
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
        lines.append(f"{int(self.charge)} {self.multiplicity}")
        if self.elements is not None:
            for i, element in enumerate(self.elements):
                lines.append(f"     {element} {self.coords[i][0]: >8.5f} {self.coords[i][1]: >8.5f} {self.coords[i][2]: >8.5f} ")
        lines.append("\n")

        return lines
    
    def print(self):
        for line in self.generate_block():
            print(line)




class GaussianReader:
    def __init__(self, filename):
        """ This is a class for reading Gaussian log files and pulling out information from them.

        Parameters
        ----------
        filename : str
            The name of the Gaussian log file to read
            
        """
        self.filename = Path(filename)
        return
    
    def check_complete(self):
        """ Check if the Gaussian calculation is complete

        This function checks if the Gaussian calculation is complete. This is done by reading the Gaussian log
        file and checking for the presence of the "Normal termination" string.

        Parameters
        ----------
        None

        Returns
        -------
        bool
            True if the calculation is complete, False otherwise

        """
        if not self.filename.exists():
            return False
        with open(self.filename, 'r') as f:
            for line in f:
                if "Normal termination" in line:
                    return True
        return False
    
    def read_log(self):
        """ Read the Gaussian log file, and extract information from it.

        This is adapted from ReadGauOutput in Tim Giese's parmutils package. This reads the 
        shebang at the end of the file, which contains information about the calculation final
        results. This works by parsing based on the \\ delimiter in this section.

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

        Notes
        -----
        .. todo::
            - Add error handling for missing data
            - Check that this reads only the FINAL geometry

        """
        atn, coords = [], []
        charge, multiplicity = 0, 1
        readflag = False
        # Open the Gaussian log file
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


