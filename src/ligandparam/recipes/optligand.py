from pathlib import Path
from typing import Optional, Union, Any

from typing_extensions import override

from ligandparam.parametrization import Recipe
from ligandparam.stages import (
    StageInitialize,
    StageDisplaceMol,
    StageNormalizeCharge,
    GaussianMinimizeRESP,
    StageLazyResp,
    StageUpdate,
    StageParmChk,
    StageLeap,
    DPMinimize,
)


class SQMLigand(Recipe):
    """This is a ligand parameterization recipe that uses Gaussian for RESP fitting.

    This script is a recipe for parameterizing a ligand using the RESP method with Gaussian. This script is designed to do the main steps of the RESP fitting process, including:
    Initializing from a pdb file, assigning atom types, running gaussian optimization, and calculating RESP charges, and then generating the final mol2/lib/frcmod files
    with RESP charges.

    Attributes:
    -----------
    in_filename : Union[Path, str]
        The input file containing the ligand structure, typically in PDB format.
    cwd : Union[Path, str]
        The current working directory where the output files will be saved.
    net_charge : int
        The net charge of the ligand.
    theory : dict
        A dictionary containing the low and high theory levels for Gaussian calculations.
    leaprc : list
        A list of leaprc files to be used in the Leap stage.
    force_gaussian_rerun : bool
        A flag indicating whether to force Gaussian reruns.
    nproc : int
        The number of processors to use for Gaussian calculations.
    mem : int
        The amount of memory (in GB) to allocate for Gaussian calculations.
    gaussian_root : Optional[Union[Path, str]]
        The root directory for Gaussian, if not set, will use the environment variable.
    gauss_exedir : Optional[Union[Path, str]]
        The directory containing the Gaussian executables, if not set, will use the environment variable.
    gaussian_binary : Optional[Union[Path, str]]
        The path to the Gaussian binary, if not set, will use the environment variable.
    gaussian_scratch : Optional[Union[Path, str]]
        The directory for Gaussian scratch files, if not set, will use the environment variable.
    kwargs : dict
        Additional keyword arguments that can be passed to the stages.
    
    Parameters:
    -----------
    in_filename : Union[Path, str]
        The input file containing the ligand structure, typically in PDB format.
    cwd : Union[Path, str]
        The current working directory where the output files will be saved.
    net_charge : int
        The net charge of the ligand.
    theory : dict, optional
        A dictionary containing the low and high theory levels for Gaussian calculations.
    leaprc : list, optional
        A list of leaprc files to be used in the Leap stage. Defaults to ["leaprc.gaff2"].
    force_gaussian_rerun : bool, optional
        A flag indicating whether to force Gaussian reruns. Defaults to False.
    nproc : int, optional
        The number of processors to use for Gaussian calculations. Defaults to 1.
    mem : int, optional
        The amount of memory (in GB) to allocate for Gaussian calculations. Defaults to 1.
    gaussian_root : Optional[Union[Path, str]], optional
        The root directory for Gaussian, if not set, will use the environment variable.
    gauss_exedir : Optional[Union[Path, str]], optional
        The directory containing the Gaussian executables, if not set, will use the environment variable.
    gaussian_binary : Optional[Union[Path, str]], optional
        The path to the Gaussian binary, if not set, will use the environment variable.
    gaussian_scratch : Optional[Union[Path, str]], optional
        The directory for Gaussian scratch files, if not set, will use the environment variable.
    kwargs : dict, optional
        Additional keyword arguments that can be passed to the stages.

    Raises:
    -------
    KeyError
        If a required option is missing from the keyword arguments.
    ValueError
        If an unknown charge model is specified in the keyword arguments.
    TypeError
        If the `theory` parameter is not a dictionary with 'low' and 'high' keys.
    AttributeError
        If a required attribute is not set during initialization.
    
    Example:
    --------
    >>> from ligandparam.recipes import LazyLigand
    >>> from pathlib import Path
    >>> recipe = LazyLigand(
        in_filename=Path("ligand.pdb"),
        cwd=Path("output_directory"),
        net_charge=0,
        theory={"low": "HF/6-31G*", "high": "PBE1PBE/6-31G*"},
        leaprc=["leaprc.gaff2"],
        force_gaussian_rerun=False,
        nproc=4,
        mem=8,
        gaussian_root=None,
        gauss_exedir=None,
        gaussian_binary=None,
        gaussian_scratch=None,
        logger="stream"
    )

    """

    @override
    def __init__(self, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs):
        super().__init__(in_filename, cwd, *args, **kwargs)
        # logger will be passed manually to each stage
        kwargs.pop("logger", None)

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
        """ Sets up the stages for the LazyLigand recipe.
        
        This method initializes the stages required for the LazyLigand recipe, including:
        - Initializing the ligand from a PDB file
        - Normalizing the charge of the ligand
        - Centering the ligand
        - Running Gaussian minimization for low theory
        - Running LazyResp for low theory
        - Running Gaussian minimization for high theory
        - Running LazyResp for high theory
        - Normalizing the charge of the ligand after high theory minimization
        - Updating the names and charges of the ligand
        - Generating the final mol2 file with RESP charges
        - Generating the frcmod and lib files for the ligand
        
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
            StageInitialize(
                "Initialize",
                main_input=self.in_filename,
                cwd=self.cwd,
                out_mol2=initial_mol2,
                net_charge=self.net_charge,
                logger=self.logger,
                **self.kwargs,
            ),
            StageNormalizeCharge(
                "Normalize1",
                main_input=initial_mol2,
                cwd=self.cwd,
                net_charge=self.net_charge,
                out_mol2=initial_mol2,
                logger=self.logger,
                **self.kwargs,
            ),
            StageDisplaceMol(
                "Centering",
                main_input=initial_mol2,
                cwd=self.cwd,
                out_mol=resp_mol2_low,
                logger=self.logger,
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
                resp_theory=self.theory["low"],
                force_gaussian_rerun=self.force_gaussian_rerun,
                out_gaussian_log=hightheory_minimization_gaussian_log,
                logger=self.logger,
                minimize=False,
                **self.kwargs,
            ),
            StageLazyResp(
                "LazyRespHigh",
                main_input=hightheory_minimization_gaussian_log,
                cwd=self.cwd,
                out_mol2=resp_mol2_high,
                net_charge=self.net_charge,
                logger=self.logger,
                **self.kwargs,
            ),
            StageNormalizeCharge(
                "Normalize2",
                main_input=resp_mol2_high,
                cwd=self.cwd,
                net_charge=self.net_charge,
                out_mol2=resp_mol2,
                logger=self.logger,
                **self.kwargs,
            ),
            StageUpdate(
                "UpdateNames",
                main_input=resp_mol2,
                cwd=self.cwd,
                source_mol2=initial_mol2,
                out_mol2=final_mol2,
                net_charge=self.net_charge,
                update_names=True,
                update_types=False,
                update_resname=True,
                logger=self.logger,
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
                logger=self.logger,
                **self.kwargs,
            ),
            StageParmChk("ParmChk", main_input=nonminimized_mol2, cwd=self.cwd, out_frcmod=frcmod,
                         logger=self.logger,
                         **self.kwargs),
            StageLeap("Leap", main_input=nonminimized_mol2, cwd=self.cwd, in_frcmod=frcmod, out_lib=lib,
                      logger=self.logger, **self.kwargs),
        ]

    @override
    def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        """ Execute the LazyLigand recipe.
        
        This method executes the LazyLigand recipe, which includes running all the stages defined in the setup method.
        
        Parameters:
        ----------
        dry_run : bool, optional
            If True, the stages will not be executed, but the commands that would be run will be printed.
        nproc : Optional[int], optional
            The number of processors to use for the calculations. If None, will use the value set in the recipe.
        mem : Optional[int], optional
            The amount of memory (in GB) to allocate for the calculations. If None, will use the value set in the recipe.
        
        Returns:
        -------
        None
        """
        self.logger.info(f"Starting the LazyLigand recipe at {self.cwd}")
        super().execute(dry_run=dry_run, nproc=nproc, mem=mem)
        self.logger.info("Done with the LazyLigand recipe")
