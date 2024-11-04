from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber
from ligandparam.io.coordinates import Remove_PDB_CONECT

class StageInitialize(AbstractStage):
    def __init__(self, name, inputoptions=None) -> None:
        """ This class is used to initialize from pdb to mol2 file using Antechamber.
        
        Parameters
        ----------
        name : str
            The name of the stage
        inputoptions : dict
            The input options
            
        """
        self.name = name
        self._parse_inputoptions(inputoptions)
        self._add_required(self.pdb_filename)
        self._add_output(f"{self.base_name}.antechamber.mol2")
        
        return
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """ Appends the stage. """
        return stage

    def _execute(self, dry_run=False):
        """ This function sets up a run in antechamber to generate a mol2 file from a pdb file with bcc charges.

        This function does a few things to get ready for a parametrization. It first removes the CONECT lines from the pdb file
        and then runs antechamber to generate a mol2 file with bcc charges and the specified atom type.
        
        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would
        
        Returns
        -------
        None
        """
        Remove_PDB_CONECT(self.pdb_filename)
        ante = Antechamber()
        ante.call(i=self.pdb_filename, fi='pdb',
                  o=self.base_name+'.antechamber.mol2', fo='mol2',
                  c='bcc', nc=self.net_charge,
                  pf='y', at=self.atom_type,
                  dry_run = dry_run)
        
    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")
    
"""
class StageSmilestoPDB(AbstractStage):
     This class is used to initialize from smiles to pdb.
    
    def __init__(self, name,=None) -> None:
        pass
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        pass
    
    def _execute(self, dry_run=False):
        pass
    
    def _clean(self):
        pass
    
    """