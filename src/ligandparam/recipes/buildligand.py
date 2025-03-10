from ligandparam.parametrization import Recipe
from ligandparam.stages import *


class BuildLigand(Recipe):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.nproc = kwargs.get("nproc", 12)
        self.mem = kwargs.get("mem", "60GB")
        self.net_charge = kwargs.get("net_charge", 0)
        self.atom_type = kwargs.get("atom_type", "gaff2")
        self.leaprc = kwargs.get("leaprc", None)
        self.target_pdb = kwargs.get("target_pdb")
        self.force_gaussian_rerun = kwargs.get("force_gaussian_rerun", False)

    def setup(self):
        self.stages = [
            StageInitialize("Initialize", inputoptions=self.inputoptions),
            StageNormalizeCharge("Normalize1",  
                                in_mol2=self.name+".antechamber.mol2", 
                                out_mol2=self.name+".antechamber.mol2",
                                inputoptions=self.inputoptions),
            GaussianMinimize("Minimize", inputoptions=self.inputoptions),
            StageGaussianRotation("Rotate", 
                                  alpha=[0, 30, 60, 90, 120, 150, 180],
                                  beta=[0, 30, 60, 90],
                                  gamma=[0],
                                  inputoptions=self.inputoptions),
            StageGaussiantoMol2("GrabGaussianCharge", inputoptions=self.inputoptions),
            StageMultiRespFit("MultiRespFit", inputoptions=self.inputoptions),
            StageUpdateCharge("UpdateCharge", 
                              in_mol2=self.name+".antechamber.mol2",
                              out_mol2=self.name+".resp.mol2",
                              charge_source="multistate",
                              inputoptions=self.inputoptions),
            StageNormalizeCharge("Normalize2",  
                                in_mol2=self.name+".resp.mol2", 
                                out_mol2=self.name+".resp.mol2",
                                inputoptions=self.inputoptions),
            StageUpdate("UpdateNames", 
                                in_mol2=self.name+'.antechamber.mol2',
                                to_update=self.name+'.resp.mol2',
                                out_mol2=self.name+'.resp.mol2', 
                                update_names=True, 
                                update_types=False,
                                update_resname=True,
                                inputoptions=self.inputoptions),
            StageUpdate("UpdateTypes", 
                                in_mol2=self.name+'.log.mol2',
                                to_update=self.name+'.resp.mol2',
                                out_mol2=self.name+'.resp.mol2', 
                                update_names=False, 
                                update_types=True,
                                inputoptions=self.inputoptions),
            StageParmChk("ParmChk", inputoptions=self.inputoptions),
            StageLeap("Leap", inputoptions=self.inputoptions),
            StageBuild("BuildGas",  build_type='gas', inputoptions=self.inputoptions),
            StageBuild("BuildAq",  build_type='aq', concentration=0.14, inputoptions=self.inputoptions),
            StageBuild("BuildTarget",  build_type='target', inputoptions=self.inputoptions)
        ]