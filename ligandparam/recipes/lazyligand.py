from ligandparam.parametrization import Recipe
from ligandparam.stages import *

class LazyLigand(Recipe):

    def __init__(self, *args, **kwargs):
        """ This is a recipe for doing a parametrization of a ligand using the RESP method without multi-state fitting.

        This recipe has a default list of stages that are run, and the stages can be disable by passing a dictionary of stages to disable to the disable stages method defined
        in the Recipe class.
        
        default_stage_list = {
            "Initialize": True,
            "Normalize1": True,
            "Minimize": True,
            "LazyResp": True,
            "UpdateCharge": True,
            "Normalize2": True,
            "UpdateNames": True,
            "UpdateTypes": True,
            "ParmChk": True,
            "Leap": True,
        }
        
        """
        super().__init__(*args, **kwargs)
        return
    
    def setup(self):
        self.stages = [
            StageInitialize("Initialize", inputoptions=self.inputoptions),
            StageNormalizeCharge("Normalize1", inputoptions=self.inputoptions, 
                    orig_mol2=self.base_name+".antechamber.mol2", 
                    new_mol2=self.base_name+".antechamber.mol2"),
            StageGaussian("Minimize", inputoptions=self.inputoptions),
            StageLazyResp("LazyResp", inputoptions=self.inputoptions),
            StageNormalizeCharge("Normalize2", inputoptions=self.inputoptions, 
                    orig_mol2=self.base_name+".resp.mol2", 
                    new_mol2=self.base_name+".resp.mol2"),
            StageUpdate("UpdateNames", inputoptions=self.inputoptions,
                    orig_mol2=self.base_name+'.antechamber.mol2',
                    to_update=self.base_name+'.resp.mol2',
                    new_mol2=self.base_name+'.resp.mol2', 
                    update_names=True, 
                    update_types=False,
                    update_resname=True),
            StageParmChk("ParmChk", inputoptions=self.inputoptions),
            StageLeap("Leap", inputoptions=self.inputoptions)
        ]