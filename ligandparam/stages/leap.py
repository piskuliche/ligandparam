import warnings
from typing import Union

import MDAnalysis as mda
from pathlib import Path

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.io.leapIO import LeapWriter
from ligandparam.interfaces import Leap


class StageLeap(AbstractStage):
    """ This is class to run a basic leap calculations on the ligand. """

    def __init__(self, stage_name: str, name: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, name, cwd, *args, **kwargs)

        self.in_frcmod = Path(self.cwd, f"{self.name.stem}.frcmod")
        self.add_required(self.in_frcmod)
        self.in_resp_mol2 = Path(self.cwd, f"{self.name.stem}.resp.mol2")
        self.add_required(self.in_resp_mol2)
        self.leaprc = kwargs.get("leaprc", ["leaprc.gaff2"])
        self.out_off = getattr(kwargs, "out_off", Path(self.cwd, f"{self.name.stem}.off"))

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
            u = mda.Universe(self.in_resp_mol2)
        resname = u.residues.resnames[0]
        # Add the leap commands
        leapgen.add_line(f"loadamberparams {self.in_frcmod.name}")
        leapgen.add_line(f"{resname} = loadmol2 {self.in_resp_mol2.name}")
        leapgen.add_line(f"saveOff {resname} {self.out_off}")
        leapgen.add_line("quit")
        # Write the leap input file
        leapgen.write(self.cwd / "tleap.param.in")
        # Call the leap program
        leap = Leap(cwd=self.cwd)
        leap.call(f=self.cwd / "tleap.param.in", dry_run=dry_run)
        return

    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")
