from pathlib import Path
from typing import Optional, Union, Any

from typing_extensions import override

from ligandparam.parametrization import Recipe
from ligandparam.stages import (
    StageInitialize,
    StageNormalizeCharge,
    GaussianMinimizeRESP,
    StageGaussiantoMol2,
    StageGaussianRotation,
    StageLazyResp,
    StageMultiRespFit,
    StageUpdateCharge,
    StageUpdate,
    StageParmChk,
    StageLeap,
)


class FreeLigand(Recipe):
    """This is a class for parametrizing a ligand that is free in solution.

    This class is designed to follow what has been the York group's best practices for parametrizing ligands.
    If your ligand is weird in any way, you should use a different class.

    This class does a parametrization using Gaussian and Antechamber, using also a multi-state RESP calculation.

    The steps are:

    1. Initialize the ligand using the PDB file.
    2. Normalize the charges to preserve neutrality.
    3. Minimize the ligand using Gaussian (a) At a low level of theory (b) At a high level of theory (c) Calculate the RESP charges using Gaussian at the low level of theory.
    4. Rotate the ligand to sample grid-based errors in resp charges
    5. Add the gaussian charges to a mol2 file.
    6. Perform a multi-state RESP fit.
    7. Update the charges in the mol2 file from the multistate fit.
    8. Normalize the charges to preserve neutrality.
    9. Update the atom types in the mol2 file to match the gaussian output.
    10. Use parmchk to generate the frcmod file.
    11. Generate the lib file with leap.

    """

    @override
    def __init__(self, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs):
        super().__init__(in_filename, cwd, *args, **kwargs)

        # required options
        for opt in ("net_charge",):
            try:
                setattr(self, opt, kwargs[opt])
                del kwargs[opt]
            except KeyError:
                raise KeyError(f"Missing {opt}")
        # required options with defaults
        # TODO: defaults should be a global singleton dict
        for opt, default_val in zip(
            ("theory", "leaprc", "force_gaussian_rerun", "nproc", "mem"),
            ({"low": "HF/6-31G*", "high": "PBE1PBE/6-31G*"}, ["leaprc.gaff2"], False, 1, 1),
        ):
            try:
                setattr(self, opt, kwargs[opt])
                del kwargs[opt]
            except KeyError:
                setattr(self, opt, default_val)

        # optional options, without defaults
        for opt in ("gaussian_root", "gauss_exedir", "gaussian_binary", "gaussian_scratch"):
            setattr(self, opt, kwargs.pop(opt, None))

        self.kwargs = kwargs

    def setup(self):
        initial_mol2 = self.cwd / f"{self.label}.initial.mol2"
        lowtheory_minimization_gaussian_log = self.cwd / f"{self.label}.lowtheory.minimization.log"
        hightheory_minimization_gaussian_log = self.cwd / f"{self.label}.hightheory.minimization.log"
        resp_mol2_low = self.cwd / f"{self.label}.minimized.lowtheory.mol2"
        resp_mol2_high = self.cwd / f"{self.label}.minimized.mol2"
        rotation_label = f"{self.label}.rotation"
        rotated_mol2 = self.cwd / f"{self.label}.rotated.mol2"
        out_respfit = self.cwd / f"respfit.charges.{self.label}"
        resp_mol2 = self.cwd / f"{self.label}.resp.mol2"
        final_mol2 = self.cwd / f"final_{self.label}.mol2"
        nonminimized_mol2 = self.cwd / f"{self.label}.mol2"
        frcmod = self.cwd / f"{self.label}.frcmod"
        lib = self.cwd / f"{self.label}.lib"

        self.stages = [
            StageInitialize(
                "Initialize",
                main_input=self.in_filename,
                cwd=self.cwd,
                out_mol2=initial_mol2,
                net_charge=self.net_charge,
                **self.kwargs,
            ),
            StageNormalizeCharge(
                "Normalize1",
                main_input=initial_mol2,
                cwd=self.cwd,
                net_charge=self.net_charge,
                out_mol2=initial_mol2,
                **self.kwargs,
            ),
            GaussianMinimizeRESP(
                "MinimizeLowTheory",
                main_input=initial_mol2,
                cwd=self.cwd,
                nproc=self.nproc,
                mem=self.mem,
                gaussian_root=self.gaussian_root,
                gauss_exedir=self.gauss_exedir,
                gaussian_binary=self.gaussian_binary,
                gaussian_scratch=self.gaussian_scratch,
                net_charge=self.net_charge,
                opt_theory=self.theory["low"],
                resp_theory=self.theory["low"],
                force_gaussian_rerun=self.force_gaussian_rerun,
                out_gaussian_log=lowtheory_minimization_gaussian_log,
                **self.kwargs,
            ),
            StageLazyResp(
                "Resp",
                main_input=lowtheory_minimization_gaussian_log,
                cwd=self.cwd,
                out_mol2=resp_mol2_low,
                net_charge=self.net_charge,
                **self.kwargs,
            ),
            GaussianMinimizeRESP(
                "MinimizeHighTheory",
                main_input=resp_mol2_low,
                cwd=self.cwd,
                nproc=self.nproc,
                mem=self.mem,
                gaussian_root=self.gaussian_root,
                gauss_exedir=self.gauss_exedir,
                gaussian_binary=self.gaussian_binary,
                gaussian_scratch=self.gaussian_scratch,
                net_charge=self.net_charge,
                opt_theory=self.theory["high"],
                resp_theory=self.theory["low"],
                force_gaussian_rerun=self.force_gaussian_rerun,
                out_gaussian_log=hightheory_minimization_gaussian_log,
                **self.kwargs,
            ),
            StageGaussiantoMol2(
                "GrabGaussianCharge",
                main_input=hightheory_minimization_gaussian_log,
                cwd=self.cwd,
                nproc=self.nproc,
                mem=self.mem,
                gaussian_root=self.gaussian_root,
                gauss_exedir=self.gauss_exedir,
                gaussian_binary=self.gaussian_binary,
                gaussian_scratch=self.gaussian_scratch,
                net_charge=self.net_charge,
                theory=self.theory,
                force_gaussian_rerun=self.force_gaussian_rerun,
                template_mol2=initial_mol2,
                out_mol2=resp_mol2_high,
                **self.kwargs,
            ),
            StageGaussianRotation(
                "Rotate",
                main_input=resp_mol2_high,
                cwd=self.cwd,
                nproc=self.nproc,
                mem=self.mem,
                gaussian_root=self.gaussian_root,
                gauss_exedir=self.gauss_exedir,
                gaussian_binary=self.gaussian_binary,
                gaussian_scratch=self.gaussian_scratch,
                net_charge=self.net_charge,
                theory=self.theory,
                force_gaussian_rerun=self.force_gaussian_rerun,
                out_gaussian_label=rotation_label,
                alpha=[0, 30, 60, 90, 120, 150, 180],
                beta=[0, 30, 60, 90],
                gamma=[0],
                **self.kwargs,
            ),
            # We know that the gaussian stages work in a "gaussianCalcs" directory. Quite hacky.
            StageMultiRespFit(
                "MultiRespFit",
                main_input=resp_mol2_high,
                cwd=self.cwd / "gaussianCalcs",
                in_gaussian_label=rotation_label,
                out_respfit=out_respfit,
                net_charge=self.net_charge,
                **self.kwargs,
            ),
            StageUpdateCharge(
                "UpdateCharge",
                main_input=resp_mol2_high,
                cwd=self.cwd,
                out_mol2=resp_mol2,
                charge_column=3,
                charge_source=out_respfit,
                net_charge=self.net_charge,
                **self.kwargs,
            ),
            StageNormalizeCharge(
                "Normalize2",
                main_input=resp_mol2,
                cwd=self.cwd,
                net_charge=self.net_charge,
                out_mol2=resp_mol2,
                **self.kwargs,
            ),
            StageUpdate(
                "UpdateNames",
                main_input=resp_mol2,
                cwd=self.cwd,
                source_mol2=initial_mol2,
                out_mol2=resp_mol2,
                update_names=True,
                update_types=False,
                update_resname=True,
                net_charge=self.net_charge,
                **self.kwargs,
            ),
            StageUpdate(
                "UpdateTypes",
                main_input=resp_mol2,
                cwd=self.cwd,
                source_mol2=initial_mol2,
                out_mol2=final_mol2,
                update_names=False,
                update_types=True,
                net_charge=self.net_charge,
                **self.kwargs,
            ),
            # Create a `nonminimized_mol2` with `initial_mol2` coordinates and  `final_mol2` charges
            StageUpdate(
                "UpdateCharges",
                main_input=initial_mol2,
                cwd=self.cwd,
                source_mol2=final_mol2,
                out_mol2=nonminimized_mol2,
                update_charges=True,
                net_charge=self.net_charge,
                **self.kwargs,
            ),
            StageParmChk("ParmChk", main_input=final_mol2, cwd=self.cwd, out_frcmod=frcmod, **self.kwargs),
            StageLeap("Leap", main_input=final_mol2, cwd=self.cwd, in_frcmod=frcmod, out_lib=lib, **self.kwargs),
        ]

    @override
    def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        self.logger.info(f"Starting the FreeLigand recipe at {self.cwd}")
        super().execute(dry_run=dry_run, nproc=nproc, mem=mem)
        self.logger.info("Done with the FreeLigand recipe")
