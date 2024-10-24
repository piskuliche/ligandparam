from ligand_param.parametrization import Recipe
from ligand_param.stages import *

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
    def setup(self):
        self.stages = [
            StageInitialize("Initialize", base_cls=self),
            StageNormalizeCharge("Normalize1", base_cls=self, 
                    orig_mol2=self.base_name+".antechamber.mol2", 
                    new_mol2=self.base_name+".antechamber.mol2"),
            StageGaussian("Minimize", base_cls=self),
            StageLazyResp("LazyResp", base_cls=self),
            StageNormalizeCharge("Normalize2", base_cls=self, 
                    orig_mol2=self.base_name+".resp.mol2", 
                    new_mol2=self.base_name+".resp.mol2"),
            StageUpdate("UpdateNames", base_cls=self,
                    orig_mol2=self.base_name+'.antechamber.mol2',
                    to_update=self.base_name+'.resp.mol2',
                    new_mol2=self.base_name+'.resp.mol2', 
                    update_names=True, 
                    update_types=False,
                    update_resname=True),
            StageParmChk("ParmChk", base_cls=self),
            StageLeap("Leap", base_cls=self)
        ]