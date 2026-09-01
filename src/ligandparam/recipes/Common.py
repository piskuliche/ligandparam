"""Shared recipe stage builders (write once)."""

from __future__ import annotations

from typing import Any, List, Optional

from ligandparam.io.Orientations import N_ORIENTATIONS_SO3_N28, legacy_euler_kwargs
from ligandparam.stages import (
    StageInitialize,
    StageNormalizeCharge,
    StageDisplaceMol,
    GaussianMinimizeRESP,
    StageGaussianToMol2,
    StageGaussianRotation,
    StageLazyResp,
    StageMultiRespFit,
    StageUpdateCharge,
    StageUpdate,
    StageLeap,
    StageParmChk,
    DPMinimize,
)


def gaussian_runtime_kwargs(recipe: Any) -> dict:
    """Common Gaussian stage kwargs from a configured recipe."""
    return {
        "nproc": recipe.nproc,
        "mem": recipe.mem,
        "gaussian_root": recipe.gaussian_root,
        "gauss_exedir": recipe.gauss_exedir,
        "gaussian_binary": recipe.gaussian_binary,
        "gaussian_scratch": recipe.gaussian_scratch,
        "force_gaussian_rerun": recipe.force_gaussian_rerun,
        "logger": recipe.logger,
    }


def rotation_stage_kwargs(recipe: Any) -> dict:
    """Kwargs for StageGaussianRotation (orientation protocol + legacy euler)."""
    rotation_kwargs = {
        **recipe.kwargs,
        "orientation_protocol": recipe.orientation_protocol,
    }
    if recipe.orientation_protocol == "legacy_euler":
        rotation_kwargs.update(legacy_euler_kwargs())
    return rotation_kwargs


def dp_minimize_stage(
    *,
    recipe: Any,
    main_input,
    out_xyz,
    out_mol2,
    steps: int,
    stage_name: str = "DPMinimize",
):
    """DeepMD / MACE minimize stage; ``steps`` is the default if kwargs omit it."""
    return DPMinimize(
        stage_name,
        main_input=main_input,
        cwd=recipe.cwd,
        out_xyz=out_xyz,
        out_mol2=out_mol2,
        model=recipe.kwargs.get("model", "deepmd_model.pb"),
        ftol=recipe.kwargs.get("ftol", 0.01),
        steps=recipe.kwargs.get("steps", steps),
        logger=recipe.logger,
    )


def init_normalize_center_stages(
    *,
    recipe: Any,
    initial_mol2,
    centered_out,
) -> List:
    """Initialize -> normalize charges -> center (DisplaceMol)."""
    return [
        StageInitialize(
            "Initialize",
            main_input=recipe.in_filename,
            cwd=recipe.cwd,
            out_mol2=initial_mol2,
            net_charge=recipe.net_charge,
            logger=recipe.logger,
            **recipe.kwargs,
        ),
        StageNormalizeCharge(
            "Normalize1",
            main_input=initial_mol2,
            cwd=recipe.cwd,
            net_charge=recipe.net_charge,
            out_mol2=initial_mol2,
            logger=recipe.logger,
            **recipe.kwargs,
        ),
        StageDisplaceMol(
            "Centering",
            main_input=initial_mol2,
            cwd=recipe.cwd,
            out_mol=centered_out,
            logger=recipe.logger,
        ),
    ]


def dual_minimize_lazy_resp_stages(
    *,
    recipe: Any,
    centered_mol2,
    low_log,
    high_log,
    resp_mol2_low,
    resp_mol2_high,
    high_minimize: bool = True,
) -> List:
    """Low + high GaussianMinimizeRESP each followed by StageLazyResp."""
    gkw = gaussian_runtime_kwargs(recipe)
    return [
        GaussianMinimizeRESP(
            "MinimizeLowTheory",
            main_input=centered_mol2,
            cwd=recipe.cwd,
            net_charge=recipe.net_charge,
            opt_theory=recipe.theory["low"],
            resp_theory=recipe.theory["low"],
            out_gaussian_log=low_log,
            minimize=recipe.kwargs.get("minimize", True),
            **gkw,
            **recipe.kwargs,
        ),
        StageLazyResp(
            "LazyRespLow",
            main_input=low_log,
            cwd=recipe.cwd,
            net_charge=recipe.net_charge,
            out_mol2=resp_mol2_low,
            logger=recipe.logger,
            **recipe.kwargs,
        ),
        GaussianMinimizeRESP(
            "MinimizeHighTheory",
            main_input=resp_mol2_low,
            cwd=recipe.cwd,
            net_charge=recipe.net_charge,
            opt_theory=recipe.theory["high"],
            resp_theory=recipe.theory["low"],
            out_gaussian_log=high_log,
            **({"minimize": high_minimize} if not high_minimize else {}),
            **gkw,
            **recipe.kwargs,
        ),
        StageLazyResp(
            "LazyRespHigh",
            main_input=high_log,
            cwd=recipe.cwd,
            out_mol2=resp_mol2_high,
            net_charge=recipe.net_charge,
            logger=recipe.logger,
            **recipe.kwargs,
        ),
    ]


def free_minimize_resp_rotation_stages(
    *,
    recipe: Any,
    centered_mol2,
    initial_mol2,
    low_log,
    high_log,
    resp_mol2_low,
    resp_mol2_high,
    rotation_label: str,
) -> List:
    """FreeLigand path: low min+RESP, high min+GaussiantoMol2, rotation."""
    gkw = gaussian_runtime_kwargs(recipe)
    rotation_kwargs = rotation_stage_kwargs(recipe)

    return [
        GaussianMinimizeRESP(
            "MinimizeLowTheory",
            main_input=centered_mol2,
            cwd=recipe.cwd,
            net_charge=recipe.net_charge,
            opt_theory=recipe.theory["low"],
            resp_theory=recipe.theory["low"],
            out_gaussian_log=low_log,
            minimize=recipe.kwargs.get("minimize", True),
            **gkw,
            **recipe.kwargs,
        ),
        StageLazyResp(
            "Resp",
            main_input=low_log,
            cwd=recipe.cwd,
            out_mol2=resp_mol2_low,
            net_charge=recipe.net_charge,
            logger=recipe.logger,
            **recipe.kwargs,
        ),
        GaussianMinimizeRESP(
            "MinimizeHighTheory",
            main_input=resp_mol2_low,
            cwd=recipe.cwd,
            net_charge=recipe.net_charge,
            opt_theory=recipe.theory["high"],
            resp_theory=recipe.theory["low"],
            out_gaussian_log=high_log,
            minimize=recipe.kwargs.get("minimize", True),
            **gkw,
            **recipe.kwargs,
        ),
        StageGaussianToMol2(
            "GrabGaussianCharge",
            main_input=high_log,
            cwd=recipe.cwd,
            net_charge=recipe.net_charge,
            theory=recipe.theory,
            template_mol2=initial_mol2,
            out_mol2=resp_mol2_high,
            **gkw,
            **recipe.kwargs,
        ),
        StageGaussianRotation(
            "Rotate",
            main_input=resp_mol2_high,
            cwd=recipe.cwd,
            net_charge=recipe.net_charge,
            theory=recipe.theory,
            out_gaussian_label=rotation_label,
            **gkw,
            **rotation_kwargs,
        ),
    ]


def dp_high_resp_rotation_stages(
    *,
    recipe: Any,
    resp_mol2_low,
    initial_mol2,
    high_log,
    resp_mol2_high,
    rotation_label: str,
) -> List:
    """DPFreeLigand Gaussian block after DeepMD: high ESP (no min) + rotation."""
    gkw = gaussian_runtime_kwargs(recipe)
    rotation_kwargs = rotation_stage_kwargs(recipe)

    return [
        GaussianMinimizeRESP(
            "MinimizeHighTheory",
            main_input=resp_mol2_low,
            cwd=recipe.cwd,
            net_charge=recipe.net_charge,
            resp_theory=recipe.theory["low"],
            out_gaussian_log=high_log,
            minimize=False,
            **gkw,
            **recipe.kwargs,
        ),
        StageGaussianToMol2(
            "GrabGaussianCharge",
            main_input=high_log,
            cwd=recipe.cwd,
            net_charge=recipe.net_charge,
            theory=recipe.theory,
            template_mol2=initial_mol2,
            out_mol2=resp_mol2_high,
            **gkw,
            **recipe.kwargs,
        ),
        StageGaussianRotation(
            "Rotate",
            main_input=resp_mol2_high,
            cwd=recipe.cwd,
            net_charge=recipe.net_charge,
            theory=recipe.theory,
            out_gaussian_label=rotation_label,
            **gkw,
            **rotation_kwargs,
        ),
    ]


def high_theory_lazy_resp_stages(
    *,
    recipe: Any,
    main_input,
    high_log,
    resp_mol2_high,
    minimize: bool = False,
) -> List:
    """Single high-theory GaussianMinimizeRESP + StageLazyResp."""
    gkw = gaussian_runtime_kwargs(recipe)
    return [
        GaussianMinimizeRESP(
            "MinimizeHighTheory",
            main_input=main_input,
            cwd=recipe.cwd,
            net_charge=recipe.net_charge,
            resp_theory=recipe.theory["low"],
            out_gaussian_log=high_log,
            minimize=minimize,
            **gkw,
            **recipe.kwargs,
        ),
        StageLazyResp(
            "LazyRespHigh",
            main_input=high_log,
            cwd=recipe.cwd,
            out_mol2=resp_mol2_high,
            net_charge=recipe.net_charge,
            logger=recipe.logger,
            **recipe.kwargs,
        ),
    ]


def multi_resp_update_stages(
    *,
    recipe: Any,
    resp_mol2_high,
    rotation_label: str,
    out_respfit,
    resp_mol2,
    initial_mol2,
    final_mol2,
    update_types: bool = False,
    normalize_input=None,
    expected_gaussian_logs: Optional[int] = None,
) -> List:
    """MultiRespFit -> UpdateCharge -> Normalize2 -> UpdateNames[/Types]."""
    if expected_gaussian_logs is None:
        expected_gaussian_logs = N_ORIENTATIONS_SO3_N28
    if normalize_input is None:
        normalize_input = resp_mol2

    stages: List = [
        StageMultiRespFit(
            "MultiRespFit",
            main_input=resp_mol2_high,
            cwd=recipe.cwd / "gaussianCalcs",
            in_gaussian_label=rotation_label,
            out_respfit=out_respfit,
            net_charge=recipe.net_charge,
            expected_gaussian_logs=expected_gaussian_logs,
            logger=recipe.logger,
            **recipe.kwargs,
        ),
        StageUpdateCharge(
            "UpdateCharge",
            main_input=resp_mol2_high,
            cwd=recipe.cwd,
            out_mol2=resp_mol2,
            charge_column=3,
            charge_source=out_respfit,
            net_charge=recipe.net_charge,
            logger=recipe.logger,
            **recipe.kwargs,
        ),
        StageNormalizeCharge(
            "Normalize2",
            main_input=normalize_input,
            cwd=recipe.cwd,
            net_charge=recipe.net_charge,
            out_mol2=resp_mol2,
            logger=recipe.logger,
            **recipe.kwargs,
        ),
    ]
    if update_types:
        stages.extend(
            [
                StageUpdate(
                    "UpdateNames",
                    main_input=resp_mol2,
                    cwd=recipe.cwd,
                    source_mol2=initial_mol2,
                    out_mol2=resp_mol2,
                    update_names=True,
                    update_types=False,
                    update_resname=True,
                    net_charge=recipe.net_charge,
                    logger=recipe.logger,
                    **recipe.kwargs,
                ),
                StageUpdate(
                    "UpdateTypes",
                    main_input=resp_mol2,
                    cwd=recipe.cwd,
                    source_mol2=initial_mol2,
                    out_mol2=final_mol2,
                    update_names=False,
                    update_types=True,
                    net_charge=recipe.net_charge,
                    logger=recipe.logger,
                    **recipe.kwargs,
                ),
            ]
        )
    else:
        stages.append(
            StageUpdate(
                "UpdateNames",
                main_input=resp_mol2,
                cwd=recipe.cwd,
                source_mol2=initial_mol2,
                out_mol2=final_mol2,
                net_charge=recipe.net_charge,
                update_names=True,
                update_types=False,
                update_resname=True,
                logger=recipe.logger,
                **recipe.kwargs,
            )
        )
    return stages


def normalize_update_names_stages(
    *,
    recipe: Any,
    resp_mol2_high,
    resp_mol2,
    initial_mol2,
    final_mol2,
) -> List:
    """Normalize2 + UpdateNames (Lazy / DP / SQM tail before parmchk)."""
    return [
        StageNormalizeCharge(
            "Normalize2",
            main_input=resp_mol2_high,
            cwd=recipe.cwd,
            net_charge=recipe.net_charge,
            out_mol2=resp_mol2,
            logger=recipe.logger,
            **recipe.kwargs,
        ),
        StageUpdate(
            "UpdateNames",
            main_input=resp_mol2,
            cwd=recipe.cwd,
            source_mol2=initial_mol2,
            out_mol2=final_mol2,
            net_charge=recipe.net_charge,
            update_names=True,
            update_types=False,
            update_resname=True,
            logger=recipe.logger,
            **recipe.kwargs,
        ),
    ]


def charge_update_parmchk_leap_stages(
    *,
    recipe: Any,
    initial_mol2,
    final_mol2,
    nonminimized_mol2,
    frcmod,
    lib,
) -> List:
    """Tail used by most ligand recipes: copy charges onto initial coords, parmchk, leap."""
    return [
        StageUpdate(
            "UpdateCharges",
            main_input=initial_mol2,
            cwd=recipe.cwd,
            source_mol2=final_mol2,
            out_mol2=nonminimized_mol2,
            update_charges=True,
            net_charge=recipe.net_charge,
            logger=recipe.logger,
            **recipe.kwargs,
        ),
        StageParmChk(
            "ParmChk",
            main_input=nonminimized_mol2,
            cwd=recipe.cwd,
            out_frcmod=frcmod,
            logger=recipe.logger,
            **recipe.kwargs,
        ),
        StageLeap(
            "Leap",
            main_input=nonminimized_mol2,
            cwd=recipe.cwd,
            in_frcmod=frcmod,
            out_lib=lib,
            logger=recipe.logger,
            **recipe.kwargs,
        ),
    ]


def rotation_label_for_recipe(recipe: Any) -> str:
    """Namespace orientation ESP labels by protocol."""
    if recipe.orientation_protocol == "so3_n28":
        return f"{recipe.label}.rotation.so3_n28"
    return f"{recipe.label}.rotation"
