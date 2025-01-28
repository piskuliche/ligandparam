from ligandparam.parametrization import Recipe
from ligandparam.stages import *

class LazyLigand(Recipe):
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
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        return
    
    def setup(self):
        self.stages = [
            StageInitialize("Initialize", inputoptions=self.inputoptions),
            StageNormalizeCharge("Normalize1", inputoptions=self.inputoptions, 
                    orig_mol2=self.name+".antechamber.mol2", 
                    new_mol2=self.name+".antechamber.mol2"),
            StageGaussian("Minimize", inputoptions=self.inputoptions),
            StageLazyResp("LazyResp", inputoptions=self.inputoptions),
            StageNormalizeCharge("Normalize2", inputoptions=self.inputoptions, 
                    orig_mol2=self.name+".resp.mol2", 
                    new_mol2=self.name+".resp.mol2"),
            StageUpdate("UpdateNames", inputoptions=self.inputoptions,
                    orig_mol2=self.name+'.antechamber.mol2',
                    to_update=self.name+'.resp.mol2',
                    new_mol2=self.name+'.resp.mol2', 
                    update_names=True, 
                    update_types=False,
                    update_resname=True),
            StageParmChk("ParmChk", inputoptions=self.inputoptions),
            StageLeap("Leap", inputoptions=self.inputoptions)
        ]