from pathlib import Path
from typing import Union, Any

from typing_extensions import override

from ligandparam.parametrization import Recipe
from ligandparam.stages import *


class LazierLigand(Recipe):
    """ This is a class for parametrizing a simple ligand using  Antechamber.

    This class is designed to do a quick parametrization of a very standard ligand. If your
    ligand is weird in any way, you should use a different class. This class does a very simple
    parametrization using Antechamber. The steps are:
    1. Use Antechamber to get bcc or abcg2 charges on a mol2 file out of the starting PDB.
    2. Check missing parameters using ParmChk.
    3. Generate the Leap input files.
    """

    @override
    def __init__(self, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs):
        super().__init__(in_filename, cwd, *args, **kwargs)

        for opt, default_val in zip(("net_charge", "nproc"), (0, 1)):
            try:
                setattr(self, opt, kwargs[opt])
                del kwargs[opt]
            except KeyError:
                setattr(self, opt, default_val)

        self.kwargs = kwargs

    def setup(self):
        nonminimized_mol2 = self.cwd / f"{self.label}.mol2"
        frcmod = self.cwd / f"{self.label}.frcmod"
        lib = self.cwd / f"{self.label}.lib"
        final_mol2 = self.cwd / f"final_{self.label}.mol2"

        self.stages = [
            StageInitialize("Initialize", main_input=self.in_filename, cwd=self.cwd, out_mol2=nonminimized_mol2,
                            **self.kwargs),
            StageParmChk("ParmChk", main_input=nonminimized_mol2, cwd=self.cwd, out_frcmod=frcmod, **self.kwargs),
            StageLeap("Leap", main_input=nonminimized_mol2, cwd=self.cwd, in_frcmod=frcmod, out_lib=lib, **self.kwargs)
            # TODO: copy `nonminimized_mol2` to `final_mol2`?
        ]

    @override
    def execute(self, dry_run=False, nproc=1, mem=512) -> Any:
        self.logger.info(f"Starting the LazierLigand recipe at {self.cwd}")
        super().execute(dry_run=False, nproc=1, mem=512)
        self.logger.info("Done with the LazierLigand recipe")
