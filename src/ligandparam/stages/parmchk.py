from typing import Union
from typing_extensions import override
from pathlib import Path

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import ParmChk
from ligandparam.utils import find_word_and_get_line


class StageParmChk(AbstractStage):
    """This is class to run parmchk on the ligand."""

    @override
    def __init__(self, stage_name: str, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, in_filename, cwd, *args, **kwargs)
        self.in_mol2 = Path(in_filename)
        self.add_required(self.in_mol2)
        self.out_frcmod = Path(kwargs["out_frcmod"])
        self.net_charge = getattr(kwargs, "net_charge", 0.0)

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """Appends the stage."""
        return stage

    def _execute(self, dry_run=False):
        """Execute the parmchk calcualtion to obtain the frcmod.

        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would

        Returns
        -------
        None

        """
        parm = ParmChk(cwd=self.cwd, logger=self.logger)
        parm.call(i=self.in_mol2, f="mol2", o=self.out_frcmod, s=2, dry_run=dry_run)

        if lines := find_word_and_get_line(self.out_frcmod, "ATTN"):
            self.logger.error(f"ATTN found in {self.out_frcmod}\n{lines}")
            raise RuntimeError(f"ATTN found in {self.out_frcmod}\n{lines}")
        return

    def _clean(self):
        """Clean the files generated during the stage."""
        raise NotImplementedError("clean method not implemented")
