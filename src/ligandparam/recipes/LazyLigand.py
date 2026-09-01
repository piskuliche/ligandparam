from pathlib import Path
from typing import Union

from typing_extensions import override

from ligandparam.Parametrization import Recipe, configure_gaussian_recipe
from ligandparam.recipes.Common import (
    charge_update_parmchk_leap_stages,
    dual_minimize_lazy_resp_stages,
    init_normalize_center_stages,
    normalize_update_names_stages,
)
from ligandparam.recipes.DihedOptions import append_dihed_twist_stage


class LazyLigand(Recipe):
    """Parameterize a ligand with Gaussian minimization and single-orientation RESP.

    Faster than :class:`FreeLigand`: skips multi-orientation ESP sampling while
    still producing mol2/lib/frcmod via Gaussian RESP and Leap.

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
    dihed_correct : bool, optional
        Recorded for ALPS. ligandparam does not run ffpopt.
    dihed_model : str, optional
        High-level model for dihedral fitting. Default ``"qdpi2"``.
    dihed_delta : int, optional
        Wavefront dihedral step in degrees (CLI ``--delta``). Default ``10``.
    dihed_fragment_config : FragmentConfig or dict, optional
        Scission fragmentation settings. Default ``None``.
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
    >>> from ligandparam.recipes import LazyLigand
    >>> recipe = LazyLigand(
    ...     "ligand.pdb", "output", net_charge=0, nproc=4, mem=8, logger="stream"
    ... )
    """

    @override
    def __init__(self, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs):
        super().__init__(in_filename, cwd, *args, **kwargs)
        configure_gaussian_recipe(self, kwargs, with_dihed=True)

    def setup(self):
        """Build the ordered LazyLigand stage list on ``self.stages``.

        Stages cover initialization, centering, low- and high-level Gaussian
        minimization with RESP, charge/name updates, and Leap library generation.
        """
        initial_mol2 = self.cwd / f"{self.label}.initial.mol2"
        centered_mol2 = self.cwd / f"{self.label}.centered.mol2"
        lowtheory_minimization_gaussian_log = self.cwd / f"{self.label}.lowtheory.minimization.log"
        hightheory_minimization_gaussian_log = self.cwd / f"{self.label}.hightheory.minimization.log"
        resp_mol2_low = self.cwd / f"{self.label}.lowtheory.mol2"
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
                centered_out=centered_mol2,
            ),
            *dual_minimize_lazy_resp_stages(
                recipe=self,
                centered_mol2=centered_mol2,
                low_log=lowtheory_minimization_gaussian_log,
                high_log=hightheory_minimization_gaussian_log,
                resp_mol2_low=resp_mol2_low,
                resp_mol2_high=resp_mol2_high,
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
        append_dihed_twist_stage(
            self.stages,
            recipe=self,
            mol2=nonminimized_mol2,
            lib=lib,
            frcmod=frcmod,
        )
