from ligand_param.stages.abstractstage import AbstractStage
from ligand_param.interfaces import Antechamber

class StageInitialize(AbstractStage):
    """ This class is used to initialize from pdb to mol2 file using Antechamber.

    Parameters
    ----------
    name : str
        Name of the stage.
    base_cls : object
        Object of the base class.

    """
    def __init__(self, name, base_cls=None) -> None:
        """ Initialize the StageInitialize class. 
        
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
        """ Execute the Gaussian calculations.
        
        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would
        
        Returns
        -------
        None
        """
        ante = Antechamber()
        ante.call(i=self.base_cls.base_name+'.pdb', fi='pdb',
                  o=self.base_cls.base_name+'.antechamber.mol2', fo='mol2',
                  c='bcc', nc=self.base_cls.net_charge,
                  pf='y', at=self.base_cls.atom_type,
                  dry_run = dry_run)
        
    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")