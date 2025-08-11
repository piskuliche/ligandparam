
"""
StageParmChk module
-------------------
This module provides the StageParmChk class for running parmchk on a ligand to generate frcmod files.
"""

from typing import Optional,  Union, Any
from typing_extensions import override
from pathlib import Path

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import ParmChk
from ligandparam.utils import find_word_and_get_line


class StageParmChk(AbstractStage):
    """
    Run parmchk on the ligand to generate frcmod files.

    Parameters
    ----------
    stage_name : str
        The name of the stage.
    main_input : Union[Path, str]
        Path to the input mol2 file.
    cwd : Union[Path, str]
        Current working directory.
    out_frcmod : str
        Path to the output frcmod file.
    net_charge : float, optional
        Net charge for the molecule (default: 0.0).

    Attributes
    ----------
    in_mol2 : Path
        Path to the input mol2 file.
    out_frcmod : Path
        Path to the output frcmod file.
    net_charge : float
        Net charge for the molecule.
    """

    @override
    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        """
        Initialize the StageParmChk class.

        Parameters
        ----------
        stage_name : str
            The name of the stage.
        main_input : Union[Path, str]
            Path to the input mol2 file.
        cwd : Union[Path, str]
            Current working directory.
        *args
            Additional positional arguments.
        **kwargs
            Additional keyword arguments.
        """
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_mol2 = Path(main_input)
        self.add_required(self.in_mol2)
        self.out_frcmod = Path(kwargs["out_frcmod"])
        self.net_charge = kwargs.get("net_charge", 0.0)

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """
        Append a stage to the current stage.

        Parameters
        ----------
        stage : AbstractStage
            The stage to append.

        Returns
        -------
        AbstractStage
            The appended stage.
        """
        return stage

    def execute(self, dry_run=False, nproc: Optional[int]=None, mem: Optional[int]=None) -> Any:
        """
        Execute the parmchk calculation to obtain the frcmod file.

        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would be run.
        nproc : int, optional
            Number of processors to use.
        mem : int, optional
            Amount of memory to use (in GB).

        Returns
        -------
        None

        Raises
        ------
        RuntimeError
            If 'ATTN' is found in the output frcmod file.
        """
        parm = ParmChk(cwd=self.cwd, logger=self.logger)
        parm.call(i=self.in_mol2, f="mol2", o=self.out_frcmod, s=2, dry_run=dry_run)

        if lines := find_word_and_get_line(self.out_frcmod, "ATTN"):
            self.logger.error(f"ATTN found in {self.out_frcmod}\n{lines}")
            raise RuntimeError(f"ATTN found in {self.out_frcmod}\n{lines}")
        return

    def _clean(self):
        """
        Clean the files generated during the stage.
        """
        raise NotImplementedError("clean method not implemented")
