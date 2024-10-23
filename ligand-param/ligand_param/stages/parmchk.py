import MDAnalysis as mda
import numpy as np

from ligand_param.stages.abstractstage import AbstractStage
from ligand_param.interfaces import ParmChk

class StageParmChk(AbstractStage):
    """ This is class to run parmchk on the ligand. """
    def __init__(self, name, base_cls=None) -> None:
        """ Initialize the StageGaussian class. 
        
        Parameters
        ----------
        name : str
            The name of the stage
        base_cls : Ligand
            The base class of the ligand
        
        """
        self.name = name
        self.base_cls = base_cls
        self.add_required(f"{self.base_cls.base_name}.resp.mol2")
        return
    

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """ Appends the stage. """
        return stage


    def _execute(self, dry_run=False):
        """ Execute the parmchk calcualtion to obtain the frcmod.
        
        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would
        
        Returns
        -------
        None

        """
        print(f"Executing {self.name} with netcharge={self.base_cls.net_charge}")
        parm = ParmChk()
        parm.call(i=self.base_cls.base_name+'.resp.mol2', f="mol2",
                  o=self.base_cls.base_name+'.frcmod', 
                  s=2, dry_run = dry_run)
        return
    
    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")