import os
import pathlib

from ligand_param.io.coordinates import Coordinates
from ligand_param.driver import Driver
from ligand_param.stages.initialize import StageInitialize
from ligand_param.stages.gaussian import StageGaussian, StageGaussianRotation, StageGaussiantoMol2
from ligand_param.stages.resp import StageLazyResp, StageMultiRespFit
from ligand_param.stages.leap import StageLeap
from ligand_param.stages.parmchk import StageParmChk
from ligand_param.stages.charge import StageUpdateCharge, StageNormalizeCharge
from ligand_param.stages.typematching import StageUpdateTypes


class Parametrization(Driver):
    """This is the base class for all parametrizations, that is a sub class of the Driver class (located in Driver.py).
    
    The rough approach to using this class is to generate a new Parametrization class, and then generate self.stages as a list of stages that you want to run.


    """


    def __init__(self, pdb_file, netcharge=0.0, atom_type="gaff2", 
                 theory_low='HF/6-31G*', theory_high='PBE1PBE/6-31G*', 
                 nproc=6, mem='8GB',
                leaprc = [],
                force_gaussian_rerun=False):
        """Initialize the class with a PDB file and a net charge.

        Parameters
        ----------
        pdb_file : str, optional
            The path to a PDB file containing the ligand structure.
        netcharge : int, optional
            The net charge of the ligand.
        atom_type : str, optional
            The atom type to use for the ligand. Default is 'gaff2'.

        """
        super().__init__()
        # Inputs
        self.pdb_filename = pdb_file
        self.net_charge = netcharge
        self.atom_type = atom_type
        self.theory={"low":theory_low, 
                     "high":theory_high}
        self.force_gaussian_rerun = force_gaussian_rerun


        # Settings for Processor Environment
        self.nproc = nproc
        self.mem = mem

        self.leaprc = []
        if not leaprc:
            self.leaprc = ['leaprc.gaff2']


        # Set the base name
        self.base_name = pathlib.Path(self.pdb_filename).stem

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
        if self.force_gaussian_rerun:
            print("Forcing Gaussian calculations to rerun.")
        else:
            print("Defaulting to NOT rerunning Gaussian calculations.")



    
class LazyLigand(Parametrization):
    """ This is a class for parametrizing a simple ligand using Gaussian and Antechamber.
    
    This class is designed to do a quick parametrization of a very standard ligand. If your
    ligand is weird in any way, you should use a different class. This class does a very simple 
    parametrization using Gaussian and Antechamber. The steps are:
    1. Initialize the ligand using the PDB file.
    2. Minimize the ligand using Gaussian at a low level of theory.
    3. Minimize the ligand using Gaussian at a high level of theory.
    4. Calculate the RESP charges using Gaussian at the low level of theory.
    5. Check the parameters using ParmChk.
    6. Generate the Leap input files.
    
    """
    def setup(self):
        self.stages = [
            StageInitialize("Initialize", base_cls=self),
            StageGaussian("Minimize", base_cls=self),
            StageLazyResp("LazyResp", base_cls=self),
            StageParmChk("ParmChk", base_cls=self),
            StageLeap("Leap", base_cls=self)
        ]

class FreeLigand(Parametrization):
    """ This is a class for parametrizing a ligand that is free in solution.
    
    This class is designed to follow what has been the York group's best practices for parametrizing ligands.
    If your ligand is weird in any way, you should use a different class. 
    
    This class does a parametrization using Gaussian and Antechamber, using also a multi-state RESP calculation.
    
    The steps are:
    1. Initialize the ligand using the PDB file.
    2. Normalize the charges to preserve neutrality.
    3. Minimize the ligand using Gaussian
    (a) At a low level of theory
    (b) At a high level of theory
    (c) Calculate the RESP charges using Gaussian at the low level of theory.
    4. Rotate the ligand to sample grid-based errors in resp charges
    5. Add the gaussian charges to a mol2 file.
    6. Perform a multi-state RESP fit.
    7. Update the charges in the mol2 file from the multistate fit.
    8. Normalize the charges to preserve neutrality.
    9. Update the atom types in the mol2 file to match the gaussian output.
    10. Use parmchk to generate the frcmod file.
    11. Generate the lib file with leap.
    """
    def setup(self):
        self.stages = [
            StageInitialize("Initialize", base_cls=self),
            #StageNormalizeCharge("Normalize", base_cls=self, 
            #                    orig_mol2=self.base_name+".antechamber.mol2", 
            #                    new_mol2=self.base_name+".antechamber.mol2"),
            StageGaussian("Minimize", base_cls=self),
            StageGaussianRotation("Rotate", base_cls=self,
                                  alpha=[0, 30, 60, 90, 120, 150, 180],
                                  beta=[0, 30, 60, 90],
                                  gamma=[0]),
            StageGaussiantoMol2("GrabGaussianCharge", base_cls=self),
            StageMultiRespFit("MultiRespFit", base_cls=self),
            StageUpdateCharge("UpdateCharge", base_cls=self,
                              orig_mol2=self.base_name+".antechamber.mol2",
                              new_mol2=self.base_name+".resp.mol2",
                              charge_source="multistate"),
            StageNormalizeCharge("Normalize", base_cls=self, 
                                orig_mol2=self.base_name+".resp.mol2", 
                                new_mol2=self.base_name+".resp.mol2"),
            StageUpdateTypes("UpdateTypes", base_cls=self,
                                orig_mol2=self.base_name+'.log.mol2',
                                to_update=self.base_name+'.resp.mol2',
                                new_mol2=self.base_name+'.resp.mol2'),
            StageParmChk("ParmChk", base_cls=self),
            StageLeap("Leap", base_cls=self)
        ]

class ProteinLigand(Parametrization):
    def new():
        pass

class RNALigand(Parametrization):
    def new():
        pass

class TestParametrization(Parametrization):
    def setup(self):
        self.stages = [
            StageInitialize("Initialize", base_cls=self),
        ]

if __name__ == "__main__":

    test = LazyLigand('FBKRP.pdb', netcharge=0)
    test.setup()
    test.execute(dry_run=False)

