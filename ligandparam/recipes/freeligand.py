from ligandparam.parametrization import Recipe
from ligandparam.stages import *

class FreeLigand(Recipe):

    def __init__(self, *args, **kwargs):
        """ This is a recipe for doing a parametrization of a ligand using the RESP method with multi-state fitting.
        
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
        self.stages = [
            StageInitialize("Initialize", inputoptions=self.inputoptions),
            StageNormalizeCharge("Normalize1", inputoptions=self.inputoptions, 
                                orig_mol2=self.base_name+".antechamber.mol2", 
                                new_mol2=self.base_name+".antechamber.mol2"),
            StageGaussian("Minimize", inputoptions=self.inputoptions),
            StageGaussianRotation("Rotate", inputoptions=self.inputoptions,
                                  alpha=[0, 30, 60, 90, 120, 150, 180],
                                  beta=[0, 30, 60, 90],
                                  gamma=[0]),
            StageGaussiantoMol2("GrabGaussianCharge", inputoptions=self.inputoptions),
            StageMultiRespFit("MultiRespFit", inputoptions=self.inputoptions),
            StageUpdateCharge("UpdateCharge", inputoptions=self.inputoptions,
                              orig_mol2=self.base_name+".antechamber.mol2",
                              new_mol2=self.base_name+".resp.mol2",
                              charge_source="multistate"),
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
            StageUpdate("UpdateTypes", inputoptions=self.inputoptions,
                                orig_mol2=self.base_name+'.log.mol2',
                                to_update=self.base_name+'.resp.mol2',
                                new_mol2=self.base_name+'.resp.mol2', 
                                update_names=False, 
                                update_types=True),
            StageParmChk("ParmChk", inputoptions=self.inputoptions),
            StageLeap("Leap", inputoptions=self.inputoptions)
        ]