from pathlib import Path
from typing import Optional, Union, Any

from typing_extensions import override

from ligandparam.parametrization import Recipe
from ligandparam.stages import StageInitialize, StageParmChk, StageLeap, StageUpdate, StageNormalizeCharge


class LazierLigand(Recipe):
    """ 
    LazierLigand is a recipe for parametrizing ligands using Antechamber and LEaP.

    It is designed to be a faster and more straightforward alternative to the FreeLigand recipe, which
    uses ABCGG2 or BCC instead of RESP for charge fitting, but retains the rest of the recipe's steps.
    
    Parameters
    ----------
    in_filename : Union[Path, str]
        The input file containing the ligand structure, typically in PDB format.
    cwd : Union[Path, str]
        The current working directory where the output files will be saved.
    net_charge : int
        The net charge of the ligand. This is a required parameter.
    nproc : int, optional
        The number of processors to use for the calculations. Default is 1.
    kwargs : dict, optional
        Additional keyword arguments that can be passed to the stages, such as `logger`, `atom_type`, etc.
    
    Attributes
    ----------
    stages : list
        A list of stages that will be executed in order to parametrize the ligand.
    logger : Any
        The logger instance used for logging messages during the execution of the recipe.
    label : str      
        A label for the recipe, typically derived from the input filename.
    in_filename : Path
        The input filename as a Path object.
    cwd : Path
        The current working directory as a Path object.
    net_charge : int
        The net charge of the ligand.
    nproc : int
        The number of processors to use for the calculations.
    kwargs : dict
        Additional keyword arguments for the stages.
    
    
    
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
        for opt, default_val in zip(("nproc",), (1,)):
            try:
                setattr(self, opt, kwargs[opt])
                del kwargs[opt]
            except KeyError:
                setattr(self, opt, default_val)

        self.kwargs = kwargs

    def setup(self):
        """ Sets up the stages for the LazierLigand recipe.
        
        This method initializes the stages required for the ligand parametrization process.
        It creates the necessary file paths and configures each stage with the appropriate parameters.
        
        The stages include:
        - StageInitialize: Initializes the ligand from the input file and generates a non-minimized mol2 file.
        - StageParmChk: Checks the parameters of the non-minimized mol2 file and generates a frcmod file.
        - StageNormalizeCharge: Normalizes the charges in the non-minimized mol2 file.
        - StageUpdate: Updates the names and types in the non-minimized mol2 file based on the fixed charge mol2 file.
        - StageUpdate: Updates the charges in the non-minimized mol2 file using the fixed charge mol2 file.
        - StageParmChk: Checks the parameters of the updated non-minimized mol2 file and generates a frcmod file.
        - StageLeap: Generates the leap input files from the non-minimized mol2 file and frcmod file.
        - StageUpdate: Optionally copies the non-minimized mol2 file to a final mol2 file.
        
        """
        nonminimized_mol2 = self.cwd / f"{self.label}.mol2"
        frcmod = self.cwd / f"{self.label}.frcmod"
        lib = self.cwd / f"{self.label}.lib"
        final_mol2 = self.cwd / f"final_{self.label}.mol2"
        fixed_charge_mol2 = self.cwd / f"fixed_charge_{self.label}.mol2"

        self.stages = [
            StageInitialize("Initialize", main_input=self.in_filename, cwd=self.cwd, out_mol2=nonminimized_mol2,
                            net_charge=self.net_charge, logger=self.logger,
                            **self.kwargs),
            StageParmChk("ParmChk", main_input=nonminimized_mol2, cwd=self.cwd, out_frcmod=frcmod,
                         logger=self.logger, **self.kwargs),
            StageNormalizeCharge(
                "Normalize2",
                main_input=nonminimized_mol2,
                cwd=self.cwd,
                net_charge=self.net_charge,
                out_mol2=fixed_charge_mol2,
                logger=self.logger,
                **self.kwargs,
            ),
            StageUpdate(
                "UpdateNames",
                main_input=nonminimized_mol2,
                cwd=self.cwd,
                source_mol2=fixed_charge_mol2,
                out_mol2=final_mol2,
                net_charge=self.net_charge,
                update_names=True,
                update_types=False,
                update_resname=True,
                logger=self.logger,
                **self.kwargs,
            ),
            # Create a `nonminimized_mol2` with `initial_mol2` coordinates and  `fixed_charge_mol2` charges
            StageUpdate(
                "UpdateCharges",
                main_input=nonminimized_mol2,
                cwd=self.cwd,
                source_mol2=fixed_charge_mol2,
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
            # TODO: copy `nonminimized_mol2` to `final_mol2`?
        ]

    @override
    def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        """Execute the LazierLigand recipe.
        This method runs the stages defined in the setup method to parametrize the ligand.
        Parameters
        ----------
        dry_run : bool, optional
            If True, the stages will not be executed, but the commands that would be run will be printed.
        nproc : Optional[int], optional
            The number of processors to use for the calculations. If None, the default value from the recipe will be used.
        mem : Optional[int], optional
            The amount of memory to allocate for the calculations. If None, the default value from the recipe will be used.
        """
        self.logger.info(f"Starting the LazierLigand recipe at {self.cwd}")
        super().execute(dry_run=False, nproc=1, mem=1)
        self.logger.info("Done with the LazierLigand recipe")
