import os

from ligand_param.io.coordinates import Coordinates
from ligand_param.driver import Driver
from ligand_param.stages.initialize import StageInitialize
from ligand_param.stages.gaussian import StageGaussian, StageGaussianRotation, StageGaussiantoMol2
from ligand_param.stages.fixcharge import StageNormalizeCharges
from ligand_param.stages.resp import StageLazyResp
from ligand_param.stages.leap import StageLeap


class Parametrization(Driver):
    """This is the base class for all parametrizations, that is a sub class of the Driver class (located in Driver.py).
    
    The rough approach to using this class is to generate a new Parametrization class, and then generate self.stages as a list of stages that you want to run.


    """


    def __init__(self, pdb_file, netcharge=0.0, atom_type=None, theory_low='HF/6-31G*', theory_high='PBE1PBE/6-31G*', nproc=6, mem='8GB', leaprc = []):
        """Initialize the class with a PDB file and a net charge.

        Parameters
        ----------
        pdb_file : str, optional
            The path to a PDB file containing the ligand structure.
        netcharge : int, optional
            The net charge of the ligand.
        atom_type : str, optional

        """
        super().__init__()
        # Inputs
        self.pdb_filename = pdb_file
        self.net_charge = netcharge
        self.atom_type = atom_type
        self.theory={"low":theory_low, 
                     "high":theory_high}
        
        # Settings for Processor Environment
        self.nproc = nproc
        self.mem = mem

        self.leaprc = []
        if not leaprc:
            self.leaprc = ['leaprc.gaff2']


        # Set the base name
        self.base_name = self.pdb_filename.strip('.pdb')

        # Generate the Coordinates
        self._generate_header(nproc, mem)
        self.coord_object = self.initial_coordinates()

        # Print out information for the user.
        self.print_info()

        return
    
    def add_leaprc(self, leaprc):
        self.leaprc.append(leaprc)
        return
    
    def initial_coordinates(self):
        try:
            coord_object = Coordinates(self.pdb_filename, filetype='pdb')
            return coord_object
        except FileExistsError:
            raise FileExistsError(f"ERROR: File {self.pdb_filename} does not exist.")
        
    def _generate_header(self, nproc, mem):
        self.header = [f'%NPROC={nproc}', f'%MEM={mem}']
        return
    
    def print_info(self):
        print(f"New Parametrization for {self.pdb_filename}")
        print("*******************************")
        print("User supplied parameters:")
        print(f"Net charge: {self.net_charge}")
        print(f"Atom type: {self.atom_type}")
        print(f"Low level QM Theory: {self.theory['low']}")
        print(f"High level QM Theory: {self.theory['high']}")
        print(f"Number of processors: {self.nproc}")
        print(f"Memory: {self.mem}")
        print("Gaussian header is:")
        print(self.header)



    
class LazyLigand(Parametrization):
    def setup(self):
        self.stages = [
            #StageInitialize("Initialize", base_cls=self),
            StageGaussian("Minimize", base_cls=self),
            StageLazyResp("LazyResp", base_cls=self),
            StageParmChk("ParmChk", base_cls=self),
            StageLeap("Leap", base_cls=self)
        ]

class FreeLigand(Parametrization):
    pass

class ProteinLigand(Parametrization):
    def new():
        pass

class RNALigand(Parametrization):
    def new():
        pass
    


if __name__ == "__main__":

    test = LazyLigand('FBKRP.pdb', netcharge=0)
    test.setup()
    test.execute(dry_run=False)

