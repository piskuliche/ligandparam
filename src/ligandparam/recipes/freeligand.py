from pathlib import Path
from typing import Union

from typing_extensions import override

from ligandparam.parametrization import Recipe
from ligandparam.stages import *


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
        for opt in ("net_charge", "nproc", "mem"):
            try:
                setattr(self, opt, kwargs[opt])
                del kwargs[opt]
            except KeyError:
                raise ValueError(f"ERROR: Please provide {opt} option as a keyword argument.")
        # required options with defaults
        # TODO: defaults should be a global singleton dict
        for opt, default_val in zip(
                ("theory", "leaprc", "force_gaussian_rerun"),
                ({"low": "HF/6-31G*", "high": "PBE1PBE/6-31G*"}, ["leaprc.gaff2"], False),
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
        minimization_gaussian_log = self.cwd / f"{self.label}.minimization.log"
        minimized_mol2 = self.cwd / f"{self.label}.minimized.mol2"
        rotation_label = f"{self.label}.rotation"
        rotated_mol2 = self.cwd / f"{self.label}.rotated.mol2"
        out_respfit = self.cwd / f"respfit.charges.{self.label}"
        resp_mol2 = self.cwd / f"{self.label}.resp.mol2"
        final_mol2 = self.cwd / f"{self.label}.mol2"
        frcmod = self.cwd / f"{self.label}.frcmod"
        lib = self.cwd / f"{self.label}.lib"

        self.stages = [
            StageInitialize(
                "Initialize", in_filename=self.in_filename, cwd=self.cwd, out_mol2=initial_mol2, **self.kwargs
            ),
            StageNormalizeCharge(
                "Normalize1",
                cwd=self.cwd,
                net_charge=self.net_charge,
                **self.kwargs,
                in_filename=initial_mol2,
                out_mol2=initial_mol2,

            ),
            StageGaussian(
                "Minimize",
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
                in_filename=initial_mol2,
                out_gaussian_log=minimization_gaussian_log,
                **self.kwargs,
            ),
            StageGaussiantoMol2(
                "GrabGaussianCharge",
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
                in_filename=minimization_gaussian_log,
                template_mol2=initial_mol2,
                out_mol2=minimized_mol2,
                **self.kwargs,
            ),
            StageGaussianRotation(
                "Rotate",
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
                in_filename=minimized_mol2,
                out_gaussian_label=rotation_label,
                alpha=[0, 30, 60, 90, 120, 150, 180],
                beta=[0, 30, 60, 90],
                gamma=[0],
                **self.kwargs,
            ),
            # We know that the gaussian stages work in a "gaussianCalcs" directory. Quite hacky.
            StageMultiRespFit(
                "MultiRespFit",
                in_filename=minimized_mol2,
                in_gaussian_label=rotation_label,
                out_respfit=out_respfit,
                cwd=self.cwd / "gaussianCalcs",
                **self.kwargs,
            ),
            StageUpdateCharge(
                "UpdateCharge",
                cwd=self.cwd,
                in_filename=minimized_mol2,
                out_mol2=resp_mol2,
                charge_column=3,
                charge_source=out_respfit,
                **self.kwargs,
            ),
            StageNormalizeCharge(
                "Normalize2",
                cwd=self.cwd,
                net_charge=self.net_charge,
                in_filename=resp_mol2,
                out_mol2=resp_mol2,
                **self.kwargs,
            ),
            StageUpdate(
                "UpdateNames",
                cwd=self.cwd,
                in_filename=resp_mol2,
                source_mol2=initial_mol2,
                out_mol2=resp_mol2,
                update_names=True,
                update_types=False,
                update_resname=True,
                **self.kwargs,
            ),
            StageUpdate(
                "UpdateTypes",
                cwd=self.cwd,
                in_filename=resp_mol2,
                source_mol2=initial_mol2,
                out_mol2=final_mol2,
                update_names=False,
                update_types=True,
                **self.kwargs,
            ),
            StageParmChk(
                "ParmChk",
                in_filename=final_mol2,
                out_frcmod=frcmod,
                cwd=self.cwd,
                **self.kwargs,
            ),
            StageLeap(
                "Leap",
                in_filename=final_mol2,
                in_frcmod=frcmod,
                out_lib=lib,
                cwd=self.cwd,
                **self.kwargs,
            ),
        ]
