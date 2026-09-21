from pathlib import Path
from typing import Union

from typing_extensions import override

from ligandparam.Parametrization import Recipe, configure_gaussian_recipe
from ligandparam.recipes.Common import (
    charge_update_parmchk_leap_stages,
    high_theory_lazy_resp_stages,
    init_normalize_center_stages,
    normalize_update_names_stages,
)


class SQMLigand(Recipe):
    """Parameterize a ligand with SQM/DeepMD-assisted minimization and RESP.

    Relaxes the ligand with DeepMD/SQM-oriented stages, then fits RESP charges
    with Gaussian and writes Leap outputs.

    Parameters
    ----------
    in_filename : path-like
        Input ligand structure (typically PDB).
    cwd : path-like
        Working directory for intermediate and output files.
    net_charge : int
        Net molecular charge.
    theory : dict, optional
        Mapping with ``low`` and ``high`` Gaussian theory levels.
    leaprc : list of str, optional
        Leaprc files for the Leap stage. Default ``["leaprc.gaff2"]``.
    force_gaussian_rerun : bool, optional
        Rerun Gaussian even if output logs already exist.
    nproc, mem : int, optional
        Gaussian processor count and memory in GB.
    gaussian_root, gauss_exedir, gaussian_binary, gaussian_scratch : optional
        Gaussian environment overrides; otherwise environment variables are used.
    **kwargs
        Extra options forwarded to stages (for example ``logger``).

    Raises
    ------
    KeyError
        If ``net_charge`` is not provided.

    Examples
    --------
    >>> from ligandparam.recipes import SQMLigand
    >>> recipe = SQMLigand(
    ...     "ligand.pdb", "output", net_charge=0, nproc=4, mem=8, logger="stream"
    ... )
    """

    @override
    def __init__(self, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs):
        super().__init__(in_filename, cwd, *args, **kwargs)
        configure_gaussian_recipe(self, kwargs)

    def setup(self):
        """Build the ordered SQMLigand stage list on ``self.stages``."""
        initial_mol2 = self.cwd / f"{self.label}.initial.mol2"
        resp_mol2_low = self.cwd / f"{self.label}.lowtheory.mol2"
        hightheory_minimization_gaussian_log = self.cwd / f"{self.label}.hightheory.minimization.log"
        resp_mol2_high = self.cwd / f"{self.label}.minimized.mol2"
        resp_mol2 = self.cwd / f"{self.label}.resp.mol2"
        final_mol2 = self.cwd / f"final_{self.label}.mol2"
        nonminimized_mol2 = self.cwd / f"{self.label}.mol2"
        frcmod = self.cwd / f"{self.label}.frcmod"
        lib = self.cwd / f"{self.label}.lib"

        self.stages = [
            *init_normalize_center_stages(
                recipe=self,
                initial_mol2=initial_mol2,
                centered_out=resp_mol2_low,
            ),
            *high_theory_lazy_resp_stages(
                recipe=self,
                main_input=resp_mol2_low,
                high_log=hightheory_minimization_gaussian_log,
                resp_mol2_high=resp_mol2_high,
                minimize=False,
            ),
            *normalize_update_names_stages(
                recipe=self,
                resp_mol2_high=resp_mol2_high,
                resp_mol2=resp_mol2,
                initial_mol2=initial_mol2,
                final_mol2=final_mol2,
            ),
            *charge_update_parmchk_leap_stages(
                recipe=self,
                initial_mol2=initial_mol2,
                final_mol2=final_mol2,
                nonminimized_mol2=nonminimized_mol2,
                frcmod=frcmod,
                lib=lib,
            ),
        ]
