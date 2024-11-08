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
    def __init__(self, name, inputoptions=None) -> None:
        """ Initialize the StageLazyResp class.
        
        Parameters
        ----------
        name : str
            The name of the stage
        inputoptions : dict
            The input options for the stage
        """
        self.name = name
        self._parse_inputoptions(inputoptions)

        self._add_required(f"gaussianCalcs/{self.base_name}.log")
        self._add_output(f"{self.base_name}.resp.mol2")
        return
    

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """ Appends the stage. """
        return stage


    def _execute(self, dry_run=False):
        """ Execute antechamber to convert the gaussian output to a mol2 file. 

        This stage will convert the gaussian output file to a mol2 file using antechamber,
        and will incorporate the resp charges calculated by gaussian into the mol2 file.
        
        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would
        """
        print(f"Executing {self.name} with netcharge={self.net_charge}")
        ante = Antechamber()
        ante.call(i=f"gaussianCalcs/{self.base_name}.log", fi='gout',
                  o=self.base_name+'.resp.mol2', fo='mol2',
                  gv=0, c='resp',
                  nc=self.net_charge,
                  at=self.atom_type, dry_run = dry_run)
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
    def __init__(self, name, inputoptions=None) -> None:
        """ Initialize the StageMultiRespFit class.

        Parameters
        ----------
        name : str
            The name of the stage
     : Ligand
            The base class of the ligand
        """
        self.name = name
        self._parse_inputoptions(inputoptions)

        self._add_required(f"{self.base_name}.log.mol2")
        self._add_required(f"gaussianCalcs/{self.base_name}.log")
        self._add_output(f"respfit.out")

        return
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """ Appends the stage. """
        return stage
    
    def _execute(self, dry_run=False):
        """Execute a multi-state respfitting calculation.

        This stage will execute a multi-state resp fitting calculation based on the gaussian output files created with 
        the gaussian stage. 

        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would

        """
        comp = parmhelper.BASH( 12 )
        model = ResidueResp( comp, 1)
        model.add_state( self.base_name, self.base_name+'.log.mol2', 
                        glob.glob("gaussianCalcs/"+self.base_name+"_*.log"), 
                        qmmask="@*" )
        model.multimolecule_fit(True)
        model.perform_fit("@*",unique_residues=False)
        with open("respfit.out", "w") as f:
            model.print_resp(fh=f)

        return

    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")