import warnings
from typing import Union

import MDAnalysis as mda
from pathlib import Path

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.io.leapIO import LeapWriter
from ligandparam.interfaces import Leap
from ligandparam.log import get_logger
from ligandparam.utils import find_word_and_get_line


class StageLeap(AbstractStage):
    """ This is class to run a basic leap calculations on the ligand. """

    def __init__(self, stage_name: str, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, in_filename, cwd, *args, **kwargs)
        self.in_mol2 = Path(in_filename)
        self.add_required(self.in_mol2)
        self.in_frcmod = kwargs["in_frcmod"]
        self.add_required(self.in_frcmod)
        self.out_lib = kwargs["out_lib"]

        self.leaprc = kwargs.get("leaprc", ["leaprc.gaff2"])

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """ Appends the stage. """
        return stage

    def _execute(self, dry_run=False):
        """ Setup and execute the leap lib file generation """
        # Generate the leap input file
        leapgen = LeapWriter("param")
        # Add the leaprc files
        for rc in self.leaprc:
            leapgen.add_leaprc(rc)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = mda.Universe(self.in_mol2)
        # Used to be: `resname = u.residues.resnames[0]`, but sometimes the resname is just a number, tripping tleap
        resname = "LIG"
        # Add the leap commands
        leapgen.add_line(f"loadamberparams {self.in_frcmod.name}")
        leapgen.add_line(f"{resname} = loadmol2 {self.in_mol2.name}")
        leapgen.add_line(f"saveOff {resname} {self.out_lib}")
        leapgen.add_line("quit")
        # Write the leap input file
        leapgen.write(self.cwd / "tleap.param.in")
        leap_log = Path(self.cwd, "leap.log")
        leap_log.unlink(missing_ok=True)
        # Call the leap program
        leap = Leap(cwd=self.cwd, logger=self.logger)
        leap.call(f=self.cwd / "tleap.param.in", dry_run=dry_run)

        if lines := find_word_and_get_line(leap_log, "Warning!"):
            self.logger.warning(f"Warning! found in {leap_log}\n{lines}")

        return

    def _clean(self):
        """Clean the files generated during the stage."""
        raise NotImplementedError("clean method not implemented")
