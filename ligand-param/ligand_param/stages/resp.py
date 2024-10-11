import MDAnalysis as mda
import numpy as np

from ligand_param.stages.abstractstage import AbstractStage
from ligand_param.interfaces import Antechamber

class StageLazyResp(AbstractStage):
    """ This class runs a 'lazy' resp calculation based on only
        a single gaussian output file. """
    def __init__(self, name, base_cls=None) -> None:
        """ Initialize the StageLazyResp class.
        
        Parameters
        ----------
        name : str
            The name of the stage
        base_cls : Ligand
            The base class of the ligand
            """
        self.name = name
        self.base_cls = base_cls
        return
    

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """ Appends the stage. """
        return stage


    def _execute(self, dry_run=False):
        """ Execute antechamber to convert the gaussian output to a mol2 file. 
        
        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would
        """
        print(f"Executing {self.name} with netcharge={self.base_cls.net_charge}")
        ante = Antechamber()
        ante.call(i=f"gaussianCalcs/{self.base_cls.base_name}.log", fi='gout',
                  o=self.base_cls.base_name+'.resp.mol2', fo='mol2',
                  gv=0, c='resp',
                  nc=self.base_cls.net_charge,
                  at=self.base_cls.atom_type, dry_run = dry_run)
        return
    
    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")
