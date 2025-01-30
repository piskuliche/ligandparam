from typing import Union

from pathlib import Path

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber
from ligandparam.io.coordinates import Remove_PDB_CONECT
from ligandparam.log import get_logger


class StageInitialize(AbstractStage):
    def __init__(self, stage_name: str, in_filename: Union[Path, str],
                 cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, in_filename, cwd, *args, **kwargs)
        self.in_pdb = Path(in_filename)
        self.add_required(self.in_pdb)
        self.out_mol2 = Path(kwargs["out_mol2"])

        self.net_charge = kwargs.get("net_charge", 0.0)
        self.atom_type = kwargs.get("atom_type", "gaff2")
        # for opt in ("net_charge", ):
        #     try:
        #         setattr(self, opt, kwargs[opt])
        #     except KeyError:
        #         raise ValueError(f"ERROR: Please provide {opt} option as a keyword argument.")
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
        Remove_PDB_CONECT(self.in_pdb)
        ante = Antechamber(cwd=self.cwd, logger=self.logger)
        ante.call(i=self.in_pdb, fi='pdb',
                  o=self.out_mol2, fo='mol2',
                  c='bcc', nc=self.net_charge,
                  pf='y', at=self.atom_type,
                  dry_run=dry_run)

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
