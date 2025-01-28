import MDAnalysis as mda
import numpy as np

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import ParmChk

class StageParmChk(AbstractStage):
    """ This is class to run parmchk on the ligand. """
    def __init__(self, name, inputoptions=None) -> None:
        """ Initialize the StageGaussian class. 
        
        Parameters
        ----------
        name : str
            The name of the stage
        inputoptions : dict
            The input options for the stage
        
        """
        self.name = name
        self._parse_inputoptions(inputoptions)
        self.add_required(f"{self.name}.resp.mol2")
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
        print(f"Executing {self.name} with netcharge={self.net_charge}")
        parm = ParmChk()
        parm.call(i=self.name+'.resp.mol2', f="mol2",
                  o=self.name+'.frcmod', 
                  s=2, dry_run = dry_run)
        return
    
    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")