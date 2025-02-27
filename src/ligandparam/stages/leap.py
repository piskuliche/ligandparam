from typing import Union

from pathlib import Path

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.io.leapIO import LeapWriter
from ligandparam.interfaces import Leap
from ligandparam.utils import find_word_and_get_line


class StageLeap(AbstractStage):
    def __init__(self, stage_name: str, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, in_filename, cwd, *args, **kwargs)
        self.in_mol2 = Path(in_filename)
        self.add_required(self.in_mol2)
        self.in_frcmod = kwargs["in_frcmod"]
        self.add_required(self.in_frcmod)
        self.out_lib = Path(kwargs["out_lib"])
        self.molname = kwargs.get("molname", "MOL")

        self.leaprc = kwargs.get("leaprc", ["leaprc.gaff2"])

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """Appends the stage."""
        return stage

    def _execute(self, dry_run=False):
        """Setup and execute the leap lib file generation"""
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
        # Call the leap program
        leap = Leap(cwd=self.cwd, logger=self.logger)
        leap.call(f=self.cwd / "tleap.param.in", dry_run=dry_run)

        if lines := find_word_and_get_line(leap_log, "Warning!"):
            self.logger.warning(f"Warning! found in {leap_log}\n{lines}")

        return

    def _clean(self):
        """Clean the files generated during the stage."""
        raise NotImplementedError("clean method not implemented")
