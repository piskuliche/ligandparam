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
        raise NotImplementedError(
            "The BuildLigand recipe is not yet implemented. Please use LazyLigand or another recipe for now."
        )

"""
    def setup(self):
        self.stages = [
            StageInitialize("Initialize",,,
            StageNormalizeCharge("Normalize1",,,
            GaussianMinimize("Minimize", inputoptions=self.inputoptions),
            StageGaussianRotation("Rotate",,,
            StageGaussiantoMol2("GrabGaussianCharge",,,
            StageMultiRespFit("MultiRespFit",,,
            StageUpdateCharge("UpdateCharge",,,
            StageNormalizeCharge("Normalize2",,,
            StageUpdate("UpdateNames",,,
            StageUpdate("UpdateTypes",,,
            StageParmChk("ParmChk",,,
            StageLeap("Leap",,,
            StageBuild("BuildGas",  build_type='gas', inputoptions=self.inputoptions),
            StageBuild("BuildAq",  build_type='aq', concentration=0.14, inputoptions=self.inputoptions),
            StageBuild("BuildTarget",  build_type='target', inputoptions=self.inputoptions)
        ]
"""