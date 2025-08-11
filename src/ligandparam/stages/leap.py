
"""
StageLeap module
----------------
This module provides the StageLeap class for generating leap library files for molecular systems.
"""

from typing import Optional,  Union, Any

from pathlib import Path
import shutil
from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.io.leapIO import LeapWriter
from ligandparam.interfaces import Leap
from ligandparam.utils import find_word_and_get_line


class StageLeap(AbstractStage):
    """
    Generate leap library files for molecular systems using AmberTools.

    Parameters
    ----------
    stage_name : str
        The name of the stage.
    main_input : Union[Path, str]
        Path to the input mol2 file.
    cwd : Union[Path, str]
        Current working directory.
    in_frcmod : str
        Path to the input frcmod file.
    out_lib : str
        Path to the output library file.
    molname : str, optional
        Molecule name (default: 'MOL').
    leaprc : list, optional
        List of leaprc configuration files (default: ['leaprc.gaff2']).

    Attributes
    ----------
    in_mol2 : Path
        Path to the input mol2 file.
    in_frcmod : str
        Path to the input frcmod file.
    out_lib : Path
        Path to the output library file.
    molname : str
        Molecule name.
    leaprc : list
        List of leaprc configuration files.
    out_pdb : Path
        Path to the output pdb file.
    """
    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        """
        Initialize the StageLeap class.

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
        self.in_frcmod = kwargs["in_frcmod"]
        self.add_required(self.in_frcmod)
        self.out_lib = Path(kwargs["out_lib"])
        self.molname = kwargs.get("molname", "MOL")
        self.leaprc = kwargs.get("leaprc", ["leaprc.gaff2"])

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
        Setup and execute the leap library file generation.

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
        """
        # Generate the leap input file
        leapgen = LeapWriter("param")
        # Add the leaprc files
        for rc in self.leaprc:
            leapgen.add_leaprc(rc)

        # Try not to overwrite an existing pdb file
        self.out_pdb = Path(self.out_lib.parent, f"{self.out_lib.stem}.pdb")
        if self.out_pdb.is_file():
            self.out_pdb = Path(self.out_lib.parent, f"tleap_{self.out_lib.stem}.pdb")

        # Add the leap commands
        leapgen.add_line(f"loadamberparams {self.in_frcmod.name}")
        leapgen.add_line(f"{self.molname} = loadmol2 {self.in_mol2.name}")
        leapgen.add_line(f'set {self.molname}.1 name "{self.molname}"')
        leapgen.add_line(f"saveOff {self.molname} {self.out_lib}")
        leapgen.add_line(f"savePDB {self.molname} {self.out_pdb}")
        leapgen.add_line("quit")
        # Write the leap input file
        leapgen.write(self.cwd / "tleap.param.in")
        leap_log = Path(self.cwd, "leap.log")
        leap_log.unlink(missing_ok=True)
        if self.out_lib.is_file():
            shutil.move(self.out_lib, self.out_lib.with_name(f"backup_{self.out_lib.name}"))
        # Call the leap program
        leap = Leap(cwd=self.cwd, logger=self.logger)
        leap.call(f=self.cwd / "tleap.param.in", dry_run=dry_run)

        if lines := find_word_and_get_line(leap_log, "Warning!"):
            self.logger.warning(f"Warning! found in {leap_log}\n{lines}")

        return

    def _clean(self):
        """
        Clean the files generated during the stage.
        """
        raise NotImplementedError("clean method not implemented")
