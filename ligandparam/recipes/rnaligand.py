from ligandparam.parametrization import Recipe
from ligandparam.stages import *

class RNALigand(Recipe):
    def __init__(self, *args, **kwargs):
        """ This is a recipe for doing a parametrization of an RNA-based ligand using the RESP method with multi-state fitting.

        This recipe has a default list of stages that are run, and the stages can be disable by passing a dictionary of stages to disable to the disable stages method defined
        in the Recipe class.
        
        default_stage_list = {
            "Initialize": True,
            "Normalize1": True,
            "Minimize": True,
            "Rotate": True,
            "GrabGaussianCharge": True,
            "MultiRespFit": True,
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
        raise NotImplementedError("This class is not yet implemented.")
        self.stages = [
            StageInitialize("Initialize", base_cls=self),
            """
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
            """
        ]