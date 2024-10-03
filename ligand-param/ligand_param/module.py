import os

from antechamber_interface import Antechamber
from gaussian import GaussianInput, GaussianWriter, GaussianReader
from coordinates import Coordinates, rotate



class Parametrization:
    """An example class that represents a simple entity."""


    def __init__(self, pdb_file, netcharge=None, atom_type=None, theory_low='HF/6-31G*', theory_high='PBE1PBE/6-31G*', nproc=6, mem='8GB', interactive=False):
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
        self.interactive=interactive

        # Set None behavior
        if self.net_charge is None:
            self.net_charge = 0.0

        if self.atom_type is None:
            self.atom_type = 'gaff2'

        # Set the base name
        self.base_name = self.pdb_filename.strip('.pdb')

        # Things set later
        self.init_gaus_run = None

        return
    
    def initialize(self):
        """
        Initialize the ligand parameterization process. This involves
        
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
                  run=False)
        
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
        self.init_gaus_run = gau
        

    def run_gaussian(self, alpha=[0,30,60,90,120,150,180], beta=[0,30,60,90], dry_run=False):
        if self.init_gaus_run is None:
            raise ValueError("Gaussian input file not initialized.")
        
        run_apply = os.system
        if dry_run:
            run_apply = print
        
        if not os.path.exists('./gaussianCalcs'):
            os.mkdir('./gaussianCalcs')
        run_apply(f'cp {self.base_name}.com ./gaussianCalcs/')
        run_apply(f'cd ./gaussianCalcs')
        run_apply(f'{self.init_gaus_run.get_run_command()}')
        for alp in alpha:
            for bet in beta:
                #TODO: add elements and header, and make sure they are consistent between steps. Probably initialized with class
                newgau = GaussianWriter(self.base_name+f'_rot_{alp}_{bet}.com')
                newcoords = rotate(coords, alpha=alp, beta=bet)
                newgau.add_block(GaussianInput(command=f"#P {self.theory['low']} OPT(CalcFC)",
                                    initial_coordinates = newcoords,
                                    elements = elements,
                                    header=header))
                run_apply(newgau.get_run_command())

        return

    def fit_resp(self):
        pass

    def correct_charge(self):
        pass

    def correct_atom_types(self):
        pass

    def write_parameters(self):
        pass 




    
    

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

