from typing import Union
from typing_extensions import override
from pathlib import Path

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import ParmChk
from ligandparam.log import get_logger


class StageParmChk(AbstractStage):
    """ This is class to run parmchk on the ligand. """

    @override
    def __init__(self, stage_name: str, name: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, name, cwd, *args, **kwargs)
        
        self.net_charge = getattr(kwargs, 'net_charge', 0.0)
        self.inmol2 = Path(self.cwd, f"{self.name.stem}.resp.mol2")
        self.outfrcmod = Path(self.cwd, f"{self.name.stem}.frcmod")
        self.add_required(self.inmol2)


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
        parm = ParmChk(cwd=self.cwd, logger=self.logger)
        parm.call(i=self.inmol2, f="mol2", o=self.outfrcmod, s=2, dry_run = dry_run)
        return
    
    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")