from pathlib import Path
from typing import Union

from typing_extensions import override

from ligandparam.Parametrization import Recipe, configure_gaussian_recipe
from ligandparam.recipes.Common import (
    charge_update_parmchk_leap_stages,
    free_minimize_resp_rotation_stages,
    init_normalize_center_stages,
    multi_resp_update_stages,
    rotation_label_for_recipe,
)
from ligandparam.recipes.DihedOptions import append_dihed_twist_stage


class FreeLigand(Recipe):
    """Parameterize a ligand with multi-orientation Gaussian RESP fitting.

    Pipeline: initialize -> normalize -> center -> low/high Gaussian minimize ->
    multi-orientation ESP (default: ``so3_n28`` quaternion pack) -> multi-RESP
    fit -> charge/name/type updates -> ``parmchk2`` / LEaP (``.frcmod`` / ``.lib``).

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
        If True, re-run Gaussian stages even when logs already show
        ``Normal termination``. Exposed on ``lig-getparam`` as ``-O`` /
        ``--force-gaussian-rerun``. Default False (skip complete jobs;
        incomplete orientation ESP jobs are re-run individually).
    orientation_protocol : {"so3_n28", "legacy_euler"}, optional
        Multi-RESP orientation set. Default ``so3_n28`` (28 quaternion-packed
        orientations). ``legacy_euler`` restores the older alpha/beta grid.
    dihed_correct : bool, optional
        Recorded for ALPS. ligandparam does not run ffpopt; after this
        recipe, run ``lig-dihed-correct`` (ALPS) on the Amber triplet.
    dihed_model : str, optional
        High-level model for dihedral fitting. Default ``"qdpi2"``.
    dihed_maxiter : int, optional
        Fit-then-rescan iterations. Default ``2``.
    dihed_delta : int, optional
        Wavefront dihedral step in degrees (same as CLI ``--delta``).
        Default ``10``.
    dihed_fragment_config : FragmentConfig or dict, optional
        Scission fragmentation settings (``FragmentConfig`` or a dict for
        ``FragmentConfig.from_dict``). Default ``None``.
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
    >>> from ligandparam.recipes import FreeLigand
    >>> recipe = FreeLigand(
    ...     "ligand.pdb", "output", net_charge=0, nproc=4, mem=8, logger="stream"
    ... )
    """

    @override
    def __init__(self, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs):
        super().__init__(in_filename, cwd, *args, **kwargs)
        configure_gaussian_recipe(
            self, kwargs, with_orientation=True, with_dihed=True
        )

    def setup(self):
        """Build the ordered FreeLigand stage list on ``self.stages``."""
        initial_mol2 = self.cwd / f"{self.label}.initial.mol2"
        centered_mol2 = self.cwd / f"{self.label}.centered.mol2"
        lowtheory_minimization_gaussian_log = self.cwd / f"{self.label}.lowtheory.minimization.log"
        hightheory_minimization_gaussian_log = self.cwd / f"{self.label}.hightheory.minimization.log"
        resp_mol2_low = self.cwd / f"{self.label}.minimized.lowtheory.mol2"
        resp_mol2_high = self.cwd / f"{self.label}.minimized.mol2"
        rotation_label = rotation_label_for_recipe(self)
        out_respfit = self.cwd / f"respfit.charges.{self.label}"
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
            *free_minimize_resp_rotation_stages(
                recipe=self,
                centered_mol2=centered_mol2,
                initial_mol2=initial_mol2,
                low_log=lowtheory_minimization_gaussian_log,
                high_log=hightheory_minimization_gaussian_log,
                resp_mol2_low=resp_mol2_low,
                resp_mol2_high=resp_mol2_high,
                rotation_label=rotation_label,
            ),
            *multi_resp_update_stages(
                recipe=self,
                resp_mol2_high=resp_mol2_high,
                rotation_label=rotation_label,
                out_respfit=out_respfit,
                resp_mol2=resp_mol2,
                initial_mol2=initial_mol2,
                final_mol2=final_mol2,
                update_types=True,
                normalize_input=resp_mol2,
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
