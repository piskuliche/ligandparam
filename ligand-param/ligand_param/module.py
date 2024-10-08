import os

from antechamber_interface import Antechamber
from gaussianIO import GaussianInput, GaussianWriter, GaussianReader
from coordinates import Coordinates
from driver import Driver
from stageinitialize import StageInitialize
from stagegaussian import StageGaussian
from stagegausrotation import StageGaussianRotation
from stagegaussiantomol2 import StageGaussiantoMol2
from stagenormalizecharges import StageNormalizeCharges


class Parametrization(Driver):
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
        super().__init__()
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
        self.coord_object = None

        self._generate_header(nproc, mem)
        self.coord_object = self.initial_coordinates()

        return
    
    def initialize(self):
        """
        Initialize the ligand parameterization process. 

        This function will generate the initial Gaussian input file and the Antechamber input file. It will also create the GaussianCalcs 
        directory in the same directory if it doesn't exist yet. 
        
        Parameters
     H -4.13200  4.46700 -0.47700 
        ----------
        None
        
        Returns
        -------
        None
        """
        raise NotImplementedError("This function should be implemented in a subclass.")
        # Grab the Coordinates from the PDB file
        try:
            self.coord_object = Coordinates(self.pdb_filename, filetype='pdb')
        except FileExistsError:
            raise FileExistsError(f"ERROR: File {self.pdb_filename} does not exist.")

        # Generate the Antechamber input file
        ante = Antechamber()
        ante.call(i=self.pdb_filename, fi = 'pdb',
                  o=self.base_name+'.mol2', fo = 'mol2',
                  c='bcc', nc=self.net_charge, 
                  pf='y', at=self.atom_type,
                  run=False)
        
        # Generate the Gaussian input file, this does a three step calculation,
        # 1. Optimize the geometry with the low level theory
        # 2. Optimize the geometry with the high level theory
        # 3. Generate the ESP charges with the low level theory
        self.header.append('%CHK='+self.base_name+'.antechamber.chk')
        gau = GaussianWriter(self.base_name+'.com')
        gau.add_block(GaussianInput(command=f"#P {self.theory['low']} OPT(CalcFC)",
                                    initial_coordinates = self.coord_object.get_coordinates(),
                                    elements = self.coord_object.get_elements(),
                                    header=self.header))
        gau.add_block(GaussianInput(command=f"#P {self.theory['high']} OPT(CalcFC) GEOM(ALLCheck) Guess(Read)", 
                                    header=self.header))
        gau.add_block(GaussianInput(command=f"#P {self.theory['low']} GEOM(AllCheck) Guess(Read) NoSymm Pop=mk IOp(6/33=2) GFInput GFPrint", 
                                    header=self.header))
        gau.write(dry_run=False)

        self.init_gaus_run = gau
        return 
        

    def run_gaussian(self, alpha=[0,30,60,90,120,150,180], beta=[0,30,60,90], dry_run=False):
        """ This function will run a series of Gaussian calculations on the initial geometry.

        Parameters
        ----------
        alpha : list, optional
            A list of angles to rotate the molecule in the alpha direction.
        beta : list, optional
            A list of angles to rotate the molecule in the beta direction.
        dry_run : bool, optional
            If True, the commands will be printed to the screen but not executed.

        Returns
        -------
        None
        """
        raise NotImplementedError("This function should be implemented in a subclass.")
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
        _, tmp_coords, _, _ = GaussianReader(self.base_name+'.log').read_log()
        self.coord_object.update_coordinates(tmp_coords)

        # Loop over the alpha and beta angles
        for alp in alpha:
            for bet in beta:
                #TODO: add elements and header, and make sure they are consistent between steps. Probably initialized with class
                newgau = GaussianWriter('gaussianCalcs/'+self.base_name+f'_rot_{alp}_{bet}.com')
                
                newgau.add_block(GaussianInput(command=f"#P {self.theory['low']} OPT(CalcFC)",
                                    initial_coordinates = self.coord_object.rotate(alpha=alp, beta=bet),
                                    elements = self.coord_object.get_elements(),
                                    header=self.header))
                newgau.write(dry_run=False)
                run_apply(newgau.get_run_command())

        return
    
    def initial_coordinates(self):
        try:
            coord_object = Coordinates(self.pdb_filename, filetype='pdb')
            return coord_object
        except FileExistsError:
            raise FileExistsError(f"ERROR: File {self.pdb_filename} does not exist.")

    def fit_resp(self):
        pass

    def correct_charge(self):
        pass

    def correct_atom_types(self):
        pass

    def write_parameters(self):
        pass 

    def _generate_header(self, nproc, mem):
        self.header = [f'%NPROC={nproc}', f'%MEM={mem}']
        return




    
    

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
    """
    print(test.pdb_filename)
    print(test.net_charge)
    print(test.atom_type)

    test.initialize()
    test.run_gaussian(dry_run=True)
    """
    #test.add_stage(StageInitialize("Initialize", base_name="FBKRP"))
    #test.add_stage(StageGaussian("Minimize", base_cls=test))
    #test.add_stage(StageGaussianRotation("Rotation1", base_cls=test))
    #test.add_stage(StageGaussiantoMol2("GaussianToMol2", base_cls=test))
    test.add_stage(StageNormalizeCharges("NormalizeCharges", mol2file=test.base_name+"log.mol2", netcharge=0.0))
    test.execute(dry_run=False)

