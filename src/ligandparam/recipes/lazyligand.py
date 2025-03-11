from pathlib import Path
from typing import Union

from typing_extensions import override

from ligandparam.parametrization import Recipe
from ligandparam.stages import (
    StageInitialize,
    StageNormalizeCharge,
    GaussianMinimizeRESP,
    StageLazyResp,
    StageUpdate,
    StageParmChk,
    StageLeap,
)


class LazyLigand(Recipe):
    """This is a class for parametrizing a simple ligand using Gaussian and Antechamber.

    This class is designed to do a quick parametrization of a very standard ligand. If your
    ligand is weird in any way, you should use a different class. This class does a very simple
    parametrization using Gaussian and Antechamber. The steps are:
    1. Initialize the ligand using the PDB file.
    2. Minimize the ligand using Gaussian at a low level of theory.
    3. Minimize the ligand using Gaussian at a high level of theory.
    4. Calculate the RESP charges using Gaussian at the low level of theory.
    5. Check the parameters using ParmChk.
    6. Generate the Leap input files.
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
                raise ValueError(f"ERROR: Please provide {opt} option as a keyword argument.")
        # required options with defaults
        # TODO: defaults should be a global singleton dict
        for opt, default_val in zip(
            ("theory", "leaprc", "force_gaussian_rerun", "nproc", "mem"),
            ({"low": "HF/6-31G*", "high": "PBE1PBE/6-31G*"}, ["leaprc.gaff2"], False, 1, 512),
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
        resp_mol2_low = self.cwd / f"{self.label}.lowtheory.mol2"
        resp_mol2_high = self.cwd / f"{self.label}.minimized.mol2"
        resp_mol2 = self.cwd / f"{self.label}.resp.mol2"
        final_mol2 = self.cwd / f"final_{self.label}.mol2"
        nonminimized_mol2 = self.cwd / f"{self.label}.mol2"
        frcmod = self.cwd / f"{self.label}.frcmod"
        lib = self.cwd / f"{self.label}.lib"

        self.stages = [
            StageInitialize("Initialize", main_input=self.in_filename, cwd=self.cwd, out_mol2=initial_mol2,
                            **self.kwargs),
            StageNormalizeCharge("Normalize1", main_input=initial_mol2, cwd=self.cwd, net_charge=self.net_charge,
                                 out_mol2=initial_mol2, **self.kwargs),
            GaussianMinimizeRESP("MinimizeLowTheory", main_input=initial_mol2, cwd=self.cwd, nproc=self.nproc,
                                 mem=self.mem, gaussian_root=self.gaussian_root, gauss_exedir=self.gauss_exedir,
                                 gaussian_binary=self.gaussian_binary, gaussian_scratch=self.gaussian_scratch,
                                 net_charge=self.net_charge, opt_theory=self.theory["low"],
                                 resp_theory=self.theory["low"], force_gaussian_rerun=self.force_gaussian_rerun,
                                 out_gaussian_log=lowtheory_minimization_gaussian_log, **self.kwargs),
            StageLazyResp("LazyRespLow", main_input=lowtheory_minimization_gaussian_log, cwd=self.cwd,
                          out_mol2=resp_mol2_low, **self.kwargs),
            GaussianMinimizeRESP("MinimizeHighTheory", main_input=resp_mol2_low, cwd=self.cwd, nproc=self.nproc,
                                 mem=self.mem, gaussian_root=self.gaussian_root, gauss_exedir=self.gauss_exedir,
                                 gaussian_binary=self.gaussian_binary, gaussian_scratch=self.gaussian_scratch,
                                 net_charge=self.net_charge, opt_theory=self.theory["high"],
                                 resp_theory=self.theory["low"], force_gaussian_rerun=self.force_gaussian_rerun,
                                 out_gaussian_log=hightheory_minimization_gaussian_log, **self.kwargs),
            StageLazyResp("LazyRespHigh", main_input=hightheory_minimization_gaussian_log, cwd=self.cwd,
                          out_mol2=resp_mol2_high, **self.kwargs),
            StageNormalizeCharge("Normalize2", main_input=resp_mol2_high, cwd=self.cwd, net_charge=self.net_charge,
                                 out_mol2=resp_mol2, **self.kwargs),
            StageUpdate("UpdateNames", main_input=resp_mol2, cwd=self.cwd, source_mol2=initial_mol2,
                        out_mol2=final_mol2, update_names=True, update_types=False, update_resname=True, **self.kwargs),
            StageParmChk("ParmChk", main_input=final_mol2, cwd=self.cwd, out_frcmod=frcmod, **self.kwargs),
            StageLeap("Leap", main_input=final_mol2, cwd=self.cwd, in_frcmod=frcmod, out_lib=lib, **self.kwargs),
            # Create a `nonminimized_mol2` with `initial_mol2` coordinates and  `final_mol2` charges
            StageUpdate("UpdateCharges", main_input=initial_mol2, cwd=self.cwd, source_mol2=final_mol2,
                        out_mol2=nonminimized_mol2, update_charges=True, **self.kwargs),
        ]

    @override
    def execute(self, dry_run=False, nproc=1, mem=512):
        self.logger.info(f"Starting the LazyLigand recipe at {self.cwd}")
        super().execute(dry_run=dry_run, nproc=nproc, mem=mem)
        self.logger.info("Done with the LazyLigand recipe")
