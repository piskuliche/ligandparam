"""Stage package - lazy exports to avoid eager optional-dep imports."""
from __future__ import annotations

from typing import Any

_EXPORTS = {
    "AbstractStage": ".AbstractStage",
    "StageLazyResp": ".Resp",
    "StageMultiRespFit": ".Resp",
    "StageParmChk": ".ParmChk",
    "StageLeap": ".Leap",
    "StageInitialize": ".Initialize",
    "GaussianMinimizeRESP": ".Gaussian",
    "StageGaussianRotation": ".Gaussian",
    "StageGaussiantoMol2": ".Gaussian",
    "StageGaussianToMol2": ".Gaussian",
    "GaussianRESP": ".Gaussian",
    "StageUpdateCharge": ".Charge",
    "StageNormalizeCharge": ".Charge",
    "StageUpdate": ".TypeMatching",
    "SDFToPDB": ".SdfConverters",
    "SDFToPDBBatch": ".SdfConverters",
    "StageSmilesToPDB": ".SmilesToPdb",
    "LigHFix": ".LigHFix",
    "StageDisplaceMol": ".DisplaceMol",
    "PDB_Name_Fixer": ".PdbNames",
    "StagePdbNameFixer": ".PdbNames",
    "DPMinimize": ".DeepMd",
    "StageSageCreate": ".GenerateSageParams",
    "StageSageToAmber": ".GenerateSageParams",
}

__all__ = list(_EXPORTS)


def __getattr__(name: str) -> Any:
    import importlib

    if name == "StageSmilestoPDB":
        from .SmilesToPdb import StageSmilesToPDB as StageSmilestoPDB

        return StageSmilestoPDB
    mod = _EXPORTS.get(name)
    if mod is not None:
        m = importlib.import_module(mod, __name__)
        return getattr(m, name)
    # ``from . import StageUtils`` would recurse through this __getattr__.
    utils = importlib.import_module(".StageUtils", __name__)
    if name == "StageUtils":
        return utils
    if hasattr(utils, name):
        return getattr(utils, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(__all__) | set(globals()))
