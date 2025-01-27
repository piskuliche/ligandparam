import glob

import MDAnalysis as mda
import numpy as np

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber

from ligandparam.multiresp import parmhelper
from ligandparam.multiresp.residueresp import ResidueResp

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
        self.cwd = getattr(self.base_cls, "cwd", None)

        self.add_required(f"./gaussianCalcs/{self.base_cls.base_name}.log")
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

class StageMultiRespFit(AbstractStage):
    """ This class runs a multi-state resp fitting calculation, based on 
        multiple gaussian output files. 
        
        TODO: Implement this class.
        TODO: Implement the clean method.
        TODO: Implement the execute method.
        TODO: Add a check that a multistate resp fit is possible. 
        """
    def __init__(self, name, base_cls=None) -> None:
        """ Initialize the StageMultiRespFit class.

        Parameters
        ----------
        name : str
            The name of the stage
        base_cls : Ligand
            The base class of the ligand
        """
        self.name = name
        self.base_cls = base_cls

        self.add_required(f"{self.base_cls.base_name}.log.mol2")
        self.add_required(f"gaussianCalcs/{self.base_cls.base_name}.log")

        return
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """ Appends the stage. """
        return stage
    
    def _execute(self, dry_run=False):
        """Execute a multi-state respfitting calculation.

        if __name__ == "__main__":
        comp = parmutils.BASH( 12 )
        model = rf.ResidueResp( comp, 1 )


        model.add_state( "$base", "$base.log.mol2", glob.glob("gaussianCalcs/$base_*.log"), qmmask="@*" )


        model.multimolecule_fit(True)
        model.perform_fit("@*",unique_residues=False)
        #model.preserve_residue_charges_by_shifting()
        model.print_resp()


        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would

        """
        comp = parmhelper.BASH( 12 )
        model = ResidueResp( comp, 1)
        model.add_state( self.base_cls.base_name, self.base_cls.base_name+'.log.mol2', 
                        glob.glob("gaussianCalcs/"+self.base_cls.base_name+"_*.log"), 
                        qmmask="@*" )
        model.multimolecule_fit(True)
        model.perform_fit("@*",unique_residues=False)
        with open("respfit.out", "w") as f:
            model.print_resp(fh=f)

        return

    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")