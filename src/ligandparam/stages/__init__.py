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
    if name == "StageSmilestoPDB":
        from .SmilesToPdb import StageSmilesToPDB as StageSmilestoPDB

        return StageSmilestoPDB
    mod = _EXPORTS.get(name)
    if mod is None:
        # Fall back to utilsstages for misc helpers historically star-imported.
        from . import StageUtils as _utils
        if hasattr(_utils, name):
            return getattr(_utils, name)
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib
    m = importlib.import_module(mod, __name__)
    return getattr(m, name)


def __dir__() -> list[str]:
    return sorted(set(__all__) | set(globals()))
