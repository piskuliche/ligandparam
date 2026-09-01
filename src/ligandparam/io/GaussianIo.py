import os

import numpy as np
from pathlib import Path


class GaussianWriter:
    def \
            __init__(self, filename):
        """Class for writing Gaussian input files.

        The filename selected will be the name of the Gaussian input file that is written to
        disk. The class also initializes the number of links to zero, and creates an empty list
        to store the GaussianInput objects.

        Parameters
        ----------
        filename : str
            The name of the file to write.
        """
        self.filename = filename
        self.out_dir = Path(filename).parent.mkdir(exist_ok=True)
        self.nlinks = 0
        self.links = []
        return
    
    def write(self, dry_run=False):
        """Write the Gaussian input file to disk.

        If ``dry_run`` is True, print the file contents instead of writing.
        LINK1 blocks are separated by the ``--Link1--`` delimiter.

        Parameters
        ----------
        dry_run : bool, optional
            If True, print the file instead of writing it to disk.
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
        """Print the Gaussian input file contents to the screen."""
        for linkno, link in enumerate(self.links):
            if linkno == 1: print("--Link1--")
            link.print()
    
    def add_block(self, block):
        """Add a GaussianInput block to this writer.

        Parameters
        ----------
        block : GaussianInput
            The GaussianInput block to append.
        """
        if isinstance(block, GaussianInput):
            self.nlinks += 1
            self.links.append(block)
    
    def get_run_command(self, extension='.com'):
        """Return a shell command string to run this Gaussian input file.

        Parameters
        ----------
        extension : str, optional
            Extension of the Gaussian input file.

        Returns
        -------
        str
            Command to run the Gaussian input file.
        """
        if extension not in self.filename:
            raise ValueError("Extension does not match filename.")
        return f"g16 < self.filename > {self.filename.strip(extension)}.log"

class GaussianInput:
    def __init__(self, command="# HF/6-31G* OPT", elements=None, initial_coordinates=None, charge=0, multiplicity=1, header=None):
        """Initialize a Gaussian input block.

        Parameters
        ----------
        command : str, optional
            The command for the Gaussian calculation.
        elements : list, optional
            A list of atomic symbols.
        initial_coordinates : np.ndarray, optional
            Atomic coordinates.
        charge : int, optional
            The charge of the molecule.
        multiplicity : int, optional
            The multiplicity of the molecule.
        header : list, optional
            Header lines (e.g. ``%NPROC``, ``%MEM``) for the input file.
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
        """Return a string representation (currently prints and returns None)."""
        print(self)
        return
    
    def generate_block(self):
        """Generate the Gaussian input block as a list of lines.

        Returns
        -------
        list of str
            Lines that make up this Gaussian input block.
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
        """Print this Gaussian input block to the screen."""
        for line in self.generate_block():
            print(line)




class GaussianReader:
    def __init__(self, filename):
        """Read Gaussian log files and extract calculation results.

        Parameters
        ----------
        filename : str
            Path to the Gaussian log file.
        """
        self.filename = Path(filename)
        return
    
    def check_complete(self):
        """Check whether the Gaussian calculation finished normally.

        Scans a short tail of the log (Normal termination is near EOF) instead
        of reading every line of large multi-orientation ESP outputs.

        Returns
        -------
        bool
            True if ``Normal termination`` is found in the log, else False.
        """
        if not self.filename.exists():
            return False
        marker = b"Normal termination"
        try:
            with open(self.filename, "rb") as f:
                f.seek(0, 2)
                size = f.tell()
                # Normal termination is written near EOF; avoid full-file scans
                # across ~28 orientation logs.
                f.seek(max(0, size - 16384))
                return marker in f.read()
        except OSError:
            return False
    
    def read_log(self):
        """Read the Gaussian log archive and extract final geometry.

        Adapted from ReadGauOutput in Tim Giese's parmutils package. Parses the
        archive section at the end of the file using the ``\\`` delimiter.

        Returns
        -------
        atn : list
            Atomic symbols.
        coords : list
            Atomic coordinates.
        charge : int
            Molecular charge.
        multiplicity : int
            Spin multiplicity.

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
        except Exception as exc:
            raise IOError(f"Error reading log file: {exc}") from exc
        
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
