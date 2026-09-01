import glob
from typing import Optional,  Union, Any
from pathlib import Path

from ligandparam.stages.AbstractStage import AbstractStage
from ligandparam.Interfaces import Antechamber
from ligandparam.io.GaussianIo import GaussianReader

from ligandparam.multiresp import ParmHelper as parmhelper
from ligandparam.multiresp.ResidueResp import ResidueResp


class StageLazyResp(AbstractStage):
    """
    Runs a 'lazy' RESP calculation based on a single Gaussian output file.

    Parameters
    ----------
    stage_name : str
        Name of the stage.
    main_input : Union[Path, str]
        Path to the input Gaussian log file.
    cwd : Union[Path, str]
        Current working directory.
    out_mol2 : str
        Path to the output mol2 file (from kwargs).
    net_charge : float, optional
        Net charge for the molecule (default is 0.0).
    atom_type : str, optional
        Atom type for antechamber (default is 'gaff2').
    molname : str, optional
        Molecule name (from kwargs).
    """

    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        """
        Initialize the StageLazyResp.

        Parameters
        ----------
        stage_name : str
            Name of the stage.
        main_input : Union[Path, str]
            Path to the input Gaussian log file.
        cwd : Union[Path, str]
            Current working directory.
        *args
            Additional positional arguments.
        **kwargs
            Additional keyword arguments. Must include 'out_mol2'.
        """
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_gaussian_log = Path(main_input)
        self.add_required(self.in_gaussian_log)
        self.out_mol2 = Path(kwargs["out_mol2"])

        self.net_charge = kwargs.get("net_charge", 0.0)
        self.atom_type = kwargs.get("atom_type", "gaff2")

        if "molname" in kwargs:
            self.additional_args = {"rn": kwargs["molname"]}
        else:
            self.additional_args = {}

    def _run(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        """
        Execute antechamber to convert the Gaussian output to a mol2 file.

        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would be run.
        nproc : int, optional
            Number of processors to use (default is None).
        mem : int, optional
            Memory to use in MB (default is None).

        Returns
        -------
        Any
            None
        """
        ante = Antechamber(cwd=self.cwd, logger=self.logger, nproc=self.nproc)
        ante.call(
            i=self.in_gaussian_log,
            fi="gout",
            o=self.out_mol2,
            fo="mol2",
            gv=0,
            c="resp",
            nc=self.net_charge,
            at=self.atom_type,
            an="no",
            dry_run=dry_run,
            **self.additional_args,
        )
        return


class StageMultiRespFit(AbstractStage):
    """Fit RESP charges from multiple Gaussian ESP orientation logs.

    Discovers logs matching ``*{in_gaussian_label}_*.log`` under ``cwd``
    (typically ``gaussianCalcs``). Works with both ``so3_n28``
    (``..._rot_q000.log``) and ``legacy_euler`` (``..._rot_0.00_0.00_0.00.log``)
    naming.

    Parameters
    ----------
    stage_name : str
        Stage name.
    main_input : path-like
        Template mol2 providing topology for the fit.
    cwd : path-like
        Directory containing the Gaussian rotation logs.
    in_gaussian_label : str
        Shared filename prefix used by :class:`StageGaussianRotation`.
    out_respfit : path-like
        Destination for the fitted charges.
    net_charge : float, optional
        Net molecular charge.
    expected_gaussian_logs : int, optional
        If set, require exactly this many completed logs before fitting.
    """

    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_gaussian_label = kwargs["in_gaussian_label"]
        self.in_mol2 = Path(main_input)
        self.in_gaussian_dir = Path(cwd)
        self.glob_str = str(self.in_gaussian_dir / f"*{self.in_gaussian_label}_*.log")
        self.out_respfit = Path(kwargs["out_respfit"])
        self.net_charge = kwargs.get("net_charge", 0.0)
        self.expected_gaussian_logs = kwargs.get("expected_gaussian_logs")

    def _run(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        """Run multi-state RESP fitting on the discovered Gaussian logs.

        Parameters
        ----------
        dry_run : bool, optional
            Currently unused; kept for stage API consistency.
        nproc, mem : optional
            Unused; kept for stage API consistency.
        """
        # Both Euler and quaternion jobs keep the "<label>_*.log" contract.
        # Sorting makes the fit deterministic (q000, q001, ...).
        gaussian_logs = sorted(glob.glob(self.glob_str))
        if not gaussian_logs:
            raise FileNotFoundError(
                f"No Gaussian rotation logs matched {self.glob_str!r}"
            )
        if (
            self.expected_gaussian_logs is not None
            and len(gaussian_logs) != self.expected_gaussian_logs
        ):
            raise RuntimeError(
                f"Expected {self.expected_gaussian_logs} Gaussian rotation logs "
                f"for {self.in_gaussian_label!r}, found {len(gaussian_logs)}. "
                "Refusing to perform a partial or mixed-orientation RESP fit."
            )
        incomplete_logs = [
            filename
            for filename in gaussian_logs
            if not GaussianReader(filename).check_complete()
        ]
        if incomplete_logs:
            names = ", ".join(Path(filename).name for filename in incomplete_logs)
            raise RuntimeError(
                f"{len(incomplete_logs)} Gaussian rotation job(s) did not "
                f"terminate normally: {names}"
            )

        comp = parmhelper.BASH(12)
        model = ResidueResp(comp, 1)
        model.add_state(
            self.in_gaussian_label,
            str(self.in_mol2),
            gaussian_logs,
            qmmask="@*",
        )
        model.multimolecule_fit(True)
        model.perform_fit("@*", unique_residues=False)
        with open(self.out_respfit, "w") as f:
            model.print_resp(fh=f)

        return

