from ligandparam.parametrization import Recipe
from ligandparam.stages import *

class FreeLigand(Recipe):
    """ This is a class for parametrizing a ligand that is free in solution.
    
    This class is designed to follow what has been the York group's best practices for parametrizing ligands.
    If your ligand is weird in any way, you should use a different class. 
    
    This class does a parametrization using Gaussian and Antechamber, using also a multi-state RESP calculation.
    
    The steps are:

    1. Initialize the ligand using the PDB file.
    2. Normalize the charges to preserve neutrality.
    3. Minimize the ligand using Gaussian (a) At a low level of theory (b) At a high level of theory (c) Calculate the RESP charges using Gaussian at the low level of theory.
    4. Rotate the ligand to sample grid-based errors in resp charges
    5. Add the gaussian charges to a mol2 file.
    6. Perform a multi-state RESP fit.
    7. Update the charges in the mol2 file from the multistate fit.
    8. Normalize the charges to preserve neutrality.
    9. Update the atom types in the mol2 file to match the gaussian output.
    10. Use parmchk to generate the frcmod file.
    11. Generate the lib file with leap.

    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        return
    def setup(self):
        self.stages = [
            StageInitialize("Initialize", inputoptions=self.inputoptions),
            StageNormalizeCharge("Normalize1", inputoptions=self.inputoptions, 
                                in_mol2=self.name+".antechamber.mol2", 
                                out_mol2=self.name+".antechamber.mol2"),
            StageGaussian("Minimize", inputoptions=self.inputoptions),
            StageGaussianRotation("Rotate", inputoptions=self.inputoptions,
                                  alpha=[0, 30, 60, 90, 120, 150, 180],
                                  beta=[0, 30, 60, 90],
                                  gamma=[0]),
            StageGaussiantoMol2("GrabGaussianCharge", inputoptions=self.inputoptions),
            StageMultiRespFit("MultiRespFit", inputoptions=self.inputoptions),
            StageUpdateCharge("UpdateCharge", inputoptions=self.inputoptions,
                              in_mol2=self.name+".antechamber.mol2",
                              out_mol2=self.name+".resp.mol2",
                              charge_source="multistate"),
            StageNormalizeCharge("Normalize2", inputoptions=self.inputoptions, 
                                in_mol2=self.name+".resp.mol2", 
                                out_mol2=self.name+".resp.mol2"),
            StageUpdate("UpdateNames", inputoptions=self.inputoptions,
                                in_mol2=self.name+'.antechamber.mol2',
                                to_update=self.name+'.resp.mol2',
                                out_mol2=self.name+'.resp.mol2', 
                                update_names=True, 
                                update_types=False,
                                update_resname=True),
            StageUpdate("UpdateTypes", inputoptions=self.inputoptions,
                                in_mol2=self.name+'.log.mol2',
                                to_update=self.name+'.resp.mol2',
                                out_mol2=self.name+'.resp.mol2', 
                                update_names=False, 
                                update_types=True),
            StageParmChk("ParmChk", inputoptions=self.inputoptions),
            StageLeap("Leap", inputoptions=self.inputoptions)
        ]