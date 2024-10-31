from ligandparam.parametrization import Recipe
from ligandparam.stages import *

class RNALigand(Recipe):
    def __init__(self, *args, **kwargs):
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