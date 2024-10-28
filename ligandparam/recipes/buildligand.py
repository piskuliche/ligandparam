from ligandparam.parametrization import Recipe
from ligandparam.stages import *


class BuildLigand(Recipe):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        return
    def setup(self):
        self.stages = [
            StageInitialize("Initialize", base_cls=self),
            StageNormalizeCharge("Normalize1", base_cls=self, 
                                orig_mol2=self.base_name+".antechamber.mol2", 
                                new_mol2=self.base_name+".antechamber.mol2"),
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
            StageUpdate("UpdateTypes", base_cls=self,
                                orig_mol2=self.base_name+'.log.mol2',
                                to_update=self.base_name+'.resp.mol2',
                                new_mol2=self.base_name+'.resp.mol2', 
                                update_names=False, 
                                update_types=True),
            StageParmChk("ParmChk", base_cls=self),
            StageLeap("Leap", base_cls=self),
            StageBuild("BuildGas", base_cls=self, build_type='gas'),
            StageBuild("BuildAq", base_cls=self, build_type='aq', concentration=0.14),
            StageBuild("BuildTarget", base_cls=self, build_type='target', target_pdb=self.target_pdb)
        ]