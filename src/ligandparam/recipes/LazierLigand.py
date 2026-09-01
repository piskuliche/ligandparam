from pathlib import Path
from typing import Union

from typing_extensions import override

from ligandparam.Parametrization import Recipe
from ligandparam.recipes.Common import charge_update_parmchk_leap_stages
from ligandparam.stages import StageInitialize, StageUpdate, StageNormalizeCharge


class LazierLigand(Recipe):
    """Fast ligand parameterization with Antechamber charges and Leap.

    Uses BCC/AM1-BCC-style charge assignment (via Antechamber) instead of
    Gaussian RESP. Suitable when speed matters more than RESP accuracy.

    Parameters
    ----------
    in_filename : path-like
        Input ligand structure (typically PDB).
    cwd : path-like
        Working directory for intermediate and output files.
    net_charge : int
        Net molecular charge.
    nproc : int, optional
        Processor count forwarded to stages. Default is 1.
    **kwargs
        Extra options forwarded to stages (for example ``logger``, ``atom_type``).

    Raises
    ------
    KeyError
        If ``net_charge`` is not provided.
    """

    @override
    def __init__(self, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs):
        super().__init__(in_filename, cwd, *args, **kwargs)
        kwargs.pop("logger", None)

        try:
            self.net_charge = kwargs.pop("net_charge")
        except KeyError as exc:
            raise KeyError("Missing net_charge") from exc
        self.nproc = kwargs.pop("nproc", 1)
        self.kwargs = kwargs

    def setup(self):
        """Build the ordered LazierLigand stage list on ``self.stages``.

        Stages cover Antechamber initialization, charge normalization, name
        updates, parmchk, and Leap library generation.
        """
        nonminimized_mol2 = self.cwd / f"{self.label}.mol2"
        frcmod = self.cwd / f"{self.label}.frcmod"
        lib = self.cwd / f"{self.label}.lib"
        final_mol2 = self.cwd / f"final_{self.label}.mol2"
        fixed_charge_mol2 = self.cwd / f"fixed_charge_{self.label}.mol2"

        self.stages = [
            StageInitialize(
                "Initialize",
                main_input=self.in_filename,
                cwd=self.cwd,
                out_mol2=nonminimized_mol2,
                net_charge=self.net_charge,
                logger=self.logger,
                **self.kwargs,
            ),
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
            *charge_update_parmchk_leap_stages(
                recipe=self,
                initial_mol2=nonminimized_mol2,
                final_mol2=fixed_charge_mol2,
                nonminimized_mol2=nonminimized_mol2,
                frcmod=frcmod,
                lib=lib,
            ),
        ]
