from pathlib import Path
from typing import Optional, Union, Any

from typing_extensions import override

from ligandparam.parametrization import Recipe
from ligandparam.stages import StageInitialize, StageParmChk, StageLeap, StageUpdate, StageNormalizeCharge


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
        # logger will be passed manually to each stage
        kwargs.pop("logger", None)

        # required options
        for opt in ("net_charge",):
            try:
                setattr(self, opt, kwargs[opt])
                del kwargs[opt]
            except KeyError:
                raise KeyError(f"Missing {opt}")
        # required options with defaults
        for opt, default_val in zip(("nproc",), (1,)):
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
        fixed_charge_mol2 = self.cwd / f"fixed_charge_{self.label}.mol2"

        self.stages = [
            StageInitialize("Initialize", main_input=self.in_filename, cwd=self.cwd, out_mol2=nonminimized_mol2,
                            net_charge=self.net_charge, logger=self.logger,
                            **self.kwargs),
            StageParmChk("ParmChk", main_input=nonminimized_mol2, cwd=self.cwd, out_frcmod=frcmod,
                         logger=self.logger, **self.kwargs),
            StageNormalizeCharge(
                "Normalize2",
                main_input=nonminimized_mol2,
                cwd=self.cwd,
                net_charge=self.net_charge,
                out_mol2=fixed_charge_mol2,
                logger=self.logger,
                **self.kwargs,
            ),
            StageUpdate(
                "UpdateNames",
                main_input=nonminimized_mol2,
                cwd=self.cwd,
                source_mol2=fixed_charge_mol2,
                out_mol2=final_mol2,
                net_charge=self.net_charge,
                update_names=True,
                update_types=False,
                update_resname=True,
                logger=self.logger,
                **self.kwargs,
            ),
            # Create a `nonminimized_mol2` with `initial_mol2` coordinates and  `fixed_charge_mol2` charges
            StageUpdate(
                "UpdateCharges",
                main_input=nonminimized_mol2,
                cwd=self.cwd,
                source_mol2=fixed_charge_mol2,
                out_mol2=nonminimized_mol2,
                update_charges=True,
                net_charge=self.net_charge,
                logger=self.logger,
                **self.kwargs,
            ),
            StageParmChk("ParmChk", main_input=nonminimized_mol2, cwd=self.cwd, out_frcmod=frcmod,
                         logger=self.logger,
                         **self.kwargs),
            StageLeap("Leap", main_input=nonminimized_mol2, cwd=self.cwd, in_frcmod=frcmod, out_lib=lib,
                      logger=self.logger, **self.kwargs),
            # TODO: copy `nonminimized_mol2` to `final_mol2`?
        ]

    @override
    def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        self.logger.info(f"Starting the LazierLigand recipe at {self.cwd}")
        super().execute(dry_run=False, nproc=1, mem=1)
        self.logger.info("Done with the LazierLigand recipe")
