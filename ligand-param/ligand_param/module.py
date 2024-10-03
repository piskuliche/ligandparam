import os

from antechamber_interface import Antechamber
from gaussian_writer import GaussianInput, GaussianWriter
from coordinates import Coordinates



class Parametrization:
    """An example class that represents a simple entity."""


    def __init__(self, pdb_file, netcharge=None, atom_type=None, theory_low='HF/6-31G*', theory_high='PBE1PBE/6-31G*', nproc=6, mem='8GB'):
        """Initialize the class with a PDB file and a net charge.

        Parameters
        ----------
        pdb_file : str, optional
            The path to a PDB file containing the ligand structure.
        netcharge : int, optional
            The net charge of the ligand.
        """
        # Inputs
        self.pdb_filename = pdb_file
        self.net_charge = netcharge
        self.atom_type = atom_type
        self.theory={"low":theory_low, 
                     "high":theory_high}
        
        self.nproc = nproc
        self.mem = mem

        # Set None behavior
        if self.net_charge is None:
            self.net_charge = 0.0

        if self.atom_type is None:
            self.atom_type = 'gaff2'

        # Set the base name
        self.base_name = self.pdb_filename.strip('.pdb')

        return
    
    def initialize(self):
        """
        Initialize the ligand parameterization process.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        coor = Coordinates(self.pdb_filename, filetype='pdb')
        coords = coor.get_coordinates()
        elements = coor.get_elements()


        ante = Antechamber()
        ante.call(i=self.pdb_filename, fi = 'pdb',
                  o=self.base_name+'.mol2', fo = 'mol2',
                  c='bcc', nc=self.net_charge, 
                  pf='y', at=self.atom_type,
                  run=True)
        
        header = [f'%NPROC={self.nproc}', f'%MEM={self.mem}', '%CHK='+self.base_name+'.antechamber.chk']
        gau = GaussianWriter(self.base_name+'.com')
        gau.add_block(GaussianInput(command=f"#P {self.theory['low']} OPT(CalcFC)",
                                    initial_coordinates = coords,
                                    elements = elements,
                                    header=header))
        gau.add_block(GaussianInput(command=f"#P {self.theory['high']} OPT(CalcFC) GEOM(ALLCheck) Guess(Read)", 
                                    header=header))
        gau.add_block(GaussianInput(command=f"#P {self.theory['low']} GEOM(AllCheck) Guess(Read) NoSymm Pop=mk IOp(6/33=2) GFInput GFPrint", 
                                    header=header))
        gau.write(dry_run=False)
    
    

class FreeLigand(Parametrization):
    pass

class ProteinLigand(Parametrization):
    def new():
        pass

class RNALigand(Parametrization):
    def new():
        pass
    


if __name__ == "__main__":

    test = FreeLigand('FBKRP.pdb', netcharge=0)

    print(test.pdb_filename)
    print(test.net_charge)
    print(test.atom_type)

    test.initialize()

