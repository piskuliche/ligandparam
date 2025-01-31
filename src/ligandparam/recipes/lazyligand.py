from pathlib import Path
from typing import Union

from typing_extensions import override

from ligandparam.parametrization import Recipe
from ligandparam.stages import *


class LazyLigand(Recipe):
    """ This is a class for parametrizing a simple ligand using Gaussian and Antechamber.
    
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
        for opt in (
            "net_charge", "nproc", "mem", "gaussian_root", "gauss_exedir", "gaussian_binary", "gaussian_scratch"):
            try:
                setattr(self, opt, kwargs[opt])
                del kwargs[opt]
            except KeyError:
                raise ValueError(f"ERROR: Please provide {opt} option as a keyword argument.")
        # required options with defaults
        # TODO: defaults should be a global singleton dict
        for opt, default_val in zip(("theory", "leaprc", "force_gaussian_rerun"),
                                ({"low": "HF/6-31G*", "high": "PBE1PBE/6-31G*"}, ["leaprc.gaff2"], False)):
            try:
                setattr(self, opt, kwargs[opt])
                del kwargs[opt]
            except KeyError:
                setattr(self, opt, default_val)

    def setup(self):
        self.stages = [
            StageInitialize("Initialize", in_filename=self.in_filename, cwd=self.cwd,
                            out_mol2=self.cwd / f"{self.label}.antechamber.mol2", **self.kwargs),
            StageNormalizeCharge("Normalize1", cwd=self.cwd, net_charge=self.net_charge, **self.kwargs,
                                 in_filename=self.cwd / f"{self.label}.antechamber.mol2",
                                 out_mol2=self.cwd / f"{self.label}.antechamber.mol2"),
            StageGaussian("Minimize", cwd=self.cwd, nproc=self.nproc, mem=self.mem,
                          gaussian_root=self.gaussian_root, gauss_exedir=self.gauss_exedir,
                          gaussian_binary=self.gaussian_binary, gaussian_scratch=self.gaussian_scratch,
                          net_charge=self.net_charge, theory=self.theory,
                          force_gaussian_rerun=self.force_gaussian_rerun,
                          in_filename=self.in_filename,
                          out_gaussian_log=self.cwd / f"{self.label}.log", **self.kwargs),
            StageLazyResp("LazyResp", cwd=self.cwd, **self.kwargs,
                          in_filename=self.cwd / f"{self.label}.log",
                          out_mol2=self.cwd / f"{self.label}.resp.mol2"),
            StageNormalizeCharge("Normalize2", cwd=self.cwd, net_charge=self.net_charge, **self.kwargs,
                                 in_filename=self.cwd / f"{self.label}.resp.mol2",
                                 out_mol2=self.cwd / f"{self.label}.resp.mol2"),
            StageUpdate("UpdateNames", cwd=self.cwd,
                        in_filename=self.cwd / f"{self.label}.antechamber.mol2",
                        to_update=self.cwd / f"{self.label}.resp.mol2",
                        out_mol2=self.cwd / f"{self.label}.resp.mol2",
                        update_names=True,
                        update_types=False,
                        update_resname=True, **self.kwargs),
            StageParmChk("ParmChk", in_filename=self.cwd / f"{self.label}.resp.mol2",
                         out_frcmod=self.cwd / f"{self.label}.frcmod", cwd=self.cwd, **self.kwargs),
            StageLeap("Leap", in_filename=self.cwd / f"{self.label}.resp.mol2",
                      in_frcmod=self.cwd / f"{self.label}.frcmod", out_lib=self.cwd / f"{self.label}.lib",
                      cwd=self.cwd, **self.kwargs)
        ]

    @override
    def execute(self, dry_run=False):
        self.logger.info(f"Starting the LazyLigand recipe at {self.cwd}")
        super().execute(dry_run=dry_run)
        self.logger.info("Done with the LazyLigand recipe")
