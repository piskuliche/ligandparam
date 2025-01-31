from pathlib import Path
from typing import Union

from typing_extensions import override

from ligandparam.parametrization import Recipe
from ligandparam.stages import *


class LazierLigand(Recipe):
    """ This is a class for parametrizing a simple ligand using  Antechamber.

    This class is designed to do a quick parametrization of a very standard ligand. If your
    ligand is weird in any way, you should use a different class. This class does a very simple
    parametrization using Gaussian and Antechamber. The steps are:
    1. Use Antechamber to get bcc charges on a mol2 file out of the starting PDB.
    2. Check the parameters using ParmChk.
    3. Generate the Leap input files.
    """

    @override
    def __init__(self, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs):
        super().__init__(in_filename, cwd, *args, **kwargs)

        # required options
        for opt in ("net_charge", "nproc"):
            try:
                setattr(self, opt, kwargs[opt])
                del kwargs[opt]
            except KeyError:
                raise ValueError(f"ERROR: Please provide {opt} option as a keyword argument.")
        self.kwargs = kwargs

    def setup(self):
        self.stages = [
            StageInitialize("Initialize", in_filename=self.in_filename, cwd=self.cwd,
                            out_mol2=self.cwd / f"{self.label}.bcc.mol2", **self.kwargs),
            StageParmChk("ParmChk", in_filename=self.cwd / f"{self.label}.bcc.mol2",
                         out_frcmod=self.cwd / f"{self.label}.frcmod", cwd=self.cwd, **self.kwargs),
            StageLeap("Leap", in_filename=self.cwd / f"{self.label}.bcc.mol2",
                      in_frcmod=self.cwd / f"{self.label}.frcmod", out_lib=self.cwd / f"{self.label}.lib",
                      cwd=self.cwd, **self.kwargs)
        ]

    @override
    def execute(self, dry_run=False):
        self.logger.info(f"Starting the LazierLigand recipe at {self.cwd}")
        super().execute(dry_run=dry_run)
        self.logger.info("Done with the LazierLigand recipe")
