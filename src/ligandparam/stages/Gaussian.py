import os
from typing import Optional,  Union, Any
import logging
import warnings
from itertools import product

import MDAnalysis as mda

from pathlib import Path
import shutil as sh

from ligandparam.stages.AbstractStage import AbstractStage
from ligandparam.io.Coordinates import Coordinates, SimpleXYZ, Mol2Writer
from ligandparam.io.GaussianIo import GaussianWriter, GaussianInput, GaussianReader
from ligandparam.io.Orientations import (
    get_quaternion_pack,
    minimum_pairwise_rotation_angle,
    quaternion_to_matrix,
)
from ligandparam.Interfaces import Gaussian, Antechamber
from ligandparam.Log import get_logger
from ligandparam.runtime.CpuBudget import split_gaussian_orientation_budget

#
logger = logging.getLogger("ligandparam.gaussian")

# Use Opt(CalcFC) below this atom count; plain Opt at or above it.
_CALCFC_MAX_ATOMS = 50

_GAUSSIAN_PATH_OPTS = (
    "gaussian_root",
    "gauss_exedir",
    "gaussian_binary",
    "gaussian_scratch",
)


def apply_gaussian_env_paths(obj, kwargs) -> None:
    """Set Gaussian installation path attributes on ``obj`` from ``kwargs``."""
    for opt in _GAUSSIAN_PATH_OPTS:
        setattr(obj, opt, kwargs.get(opt, ""))
    if getattr(obj, "gaussian_binary", None) is None:
        obj.gaussian_binary = "g16"


def _gaussian_log_is_complete(path: Path | str | None) -> bool:
    """Return True if ``path`` exists and ends with Normal termination."""
    if path is None:
        return False
    return GaussianReader(Path(path)).check_complete()


def _should_skip_gaussian_job(
    *,
    force_rerun: bool,
    final_log: Path | str | None,
    cwd_log: Path | str | None = None,
    logger: logging.Logger | None = None,
    promote_cwd_to_final: bool = True,
) -> bool:
    """Decide whether an existing Gaussian job can be skipped on resume.

    Checks ``final_log`` first, then ``cwd_log`` (e.g. ``gaussianCalcs/*.log``
    before it was moved). Incomplete logs trigger a re-run of that job only.
    Pass ``force_rerun=True`` (CLI ``-O``) to ignore complete logs.
    """
    log = logger or logging.getLogger("ligandparam.gaussian")
    final_path = Path(final_log) if final_log is not None else None
    cwd_path = Path(cwd_log) if cwd_log is not None else None

    if force_rerun:
        log.info(
            "force_gaussian_rerun (-O): ignoring existing Gaussian logs and re-running"
        )
        return False

    if final_path is not None and _gaussian_log_is_complete(final_path):
        log.info("Skipping Gaussian job (already complete): %s", final_path)
        return True

    if cwd_path is not None and _gaussian_log_is_complete(cwd_path):
        log.info("Skipping Gaussian job (already complete): %s", cwd_path)
        if (
            promote_cwd_to_final
            and final_path is not None
            and cwd_path.resolve() != final_path.resolve()
        ):
            try:
                final_path.parent.mkdir(parents=True, exist_ok=True)
                sh.copy2(cwd_path, final_path)
                log.info("Promoted complete Gaussian log -> %s", final_path)
            except OSError as exc:
                log.warning(
                    "Could not promote complete log %s -> %s (%s)",
                    cwd_path,
                    final_path,
                    exc,
                )
        return True

    for path in (final_path, cwd_path):
        if path is not None and path.exists() and not _gaussian_log_is_complete(path):
            log.info(
                "Incomplete Gaussian log found (%s); will re-run this job",
                path,
            )
    return False


def _gaussian_link0_header(nproc, mem, chk=None) -> list[str]:
    """``%NPROC`` / ``%MEM`` Link0 lines, optional ``%chk``."""
    header = [f"%NPROC={nproc}", f"%MEM={mem}GB"]
    if chk:
        header.append(f"%chk={chk}")
    return header


def _run_gaussian_com_and_promote(stage, *, dry_run: bool) -> None:
    """Run ``stage.in_com`` under ``gaussian_cwd`` and move the log out."""
    Gaussian(
        cwd=stage.gaussian_cwd,
        logger=stage.logger,
        gaussian_root=stage.gaussian_root,
        gauss_exedir=stage.gauss_exedir,
        gaussian_binary=stage.gaussian_binary,
        gaussian_scratch=stage.gaussian_scratch,
    ).call(
        inp_pipe=stage.in_com.name,
        out_pipe=stage.out_log.name,
        dry_run=dry_run,
    )
    if not dry_run:
        sh.move(stage.out_log, stage.out_gaussian_log)


def _orientation_id_from_paths(in_com: str | Path, out_log: str | Path) -> str:
    """Stable board id: ``q012`` or ``0.00_30.00_0.00`` from ``*_rot_<id>.*``."""
    stem = Path(in_com).stem
    marker = "_rot_"
    if marker in stem:
        return stem.split(marker, 1)[1]
    return Path(out_log).stem


def _run_gaussian_rotation_job(payload: dict) -> dict:
    """Run one rotation ESP job (spawn-pool worker; must be picklable)."""
    from ligandparam.runtime.ProgressBoard import JobProgressStore

    cwd = Path(payload["cwd"])
    in_com = payload["in_com"]
    out_log = payload["out_log"]
    force = bool(payload.get("force", False))
    dry_run = bool(payload.get("dry_run", False))
    log_path = cwd / out_log
    job_id = payload.get("job_id") or _orientation_id_from_paths(in_com, out_log)
    store = None
    status_path = payload.get("status_path")
    if status_path:
        store = JobProgressStore(
            status_path,
            collection_key="orientations",
            id_header="Angle",
            title="Gaussian orientation ESP - live status",
        )

    def _set(**kwargs):
        if store is not None:
            store.update(job_id, **kwargs)

    if not force and GaussianReader(log_path).check_complete():
        _set(status="skipped", stage="finished", detail=f"already complete | {out_log}")
        return {"in_com": in_com, "status": "skipped", "job_id": job_id}

    _set(
        status="running",
        stage="gaussian",
        detail=f"{out_log} | %NProc={payload.get('job_nproc', '?')} %MEM={payload.get('job_mem', '?')}GB",
        log_path=str(log_path),
    )
    try:
        gau = Gaussian(
            cwd=cwd,
            gaussian_root=payload.get("gaussian_root", ""),
            gauss_exedir=payload.get("gauss_exedir", ""),
            gaussian_binary=payload.get("gaussian_binary", "g16"),
            gaussian_scratch=payload.get("gaussian_scratch", ""),
            logger=logging.getLogger("ligandparam.gaussian.worker"),
        )
        stem = Path(in_com).stem
        gau.call(
            inp_pipe=in_com,
            out_pipe=out_log,
            dry_run=dry_run,
            script_name=f"_gau_{stem}.sh",
            scratch=str(cwd / "tmp" / f"scratch_{stem}"),
        )
        if not dry_run and not GaussianReader(log_path).check_complete():
            raise RuntimeError(f"Gaussian did not complete normally: {log_path}")
        _set(status="done", stage="finished", detail=f"ok | {out_log}")
        return {"in_com": in_com, "status": "ok", "job_id": job_id}
    except Exception as exc:
        _set(
            status="failed",
            stage="failed",
            detail=type(exc).__name__,
            error=str(exc)[:500],
        )
        raise


def _gaussian_opt_keyword(n_atoms: int) -> str:
    """Choose ``Opt(CalcFC)`` or ``Opt`` from ligand size.

    Opt=CalcFC computes the full force-constant matrix (Hessian) at the
    initial geometry. Plain Opt starts from an inexpensive approximate
    Hessian and updates it using gradients from later optimization steps.

    The additional cost of CalcFC is roughly one frequency calculation at
    the same method and basis set:

        T_Opt ~= n_steps * T_gradient
        T_Opt(CalcFC) ~= T_Hessian + n_steps' * T_gradient

    CalcFC is worthwhile only when the better starting Hessian saves enough
    optimization steps to offset T_Hessian.

    Scaling (M basis functions, N atoms):

    - The Cartesian Hessian has (3N)^2 elements -> storage scales as O(N^2).
    - For HF and many DFT methods, Gaussian has analytical second
      derivatives. Formal scaling is broadly similar to the gradient, but
      with a much larger prefactor, memory, and disk footprint.
    - Without analytical seconds, a numerical Hessian may need up to ~6N
      gradient calculations and becomes prohibitive quickly.
    - Post-HF Hessians (MP2 and higher) get expensive at much smaller
      sizes than ordinary DFT Hessians.

    Policy here: ``Opt(CalcFC)`` when ``N < 50``, otherwise ``Opt``.
    """
    if n_atoms < _CALCFC_MAX_ATOMS:
        return "Opt(CalcFC)"
    return "Opt"


class GaussianMinimizeRESP(AbstractStage):
    """
    Run a basic Gaussian calculation on the ligand, including minimization and ESP calculation for RESP charges.

    Parameters
    ----------
    stage_name : str
        The name of the stage.
    main_input : Union[Path, str]
        Path to the input mol2 file.
    cwd : Union[Path, str]
        Current working directory.
    out_gaussian_log : str
        Path to the output Gaussian log file.
    opt_theory : str, optional
        Theory for optimization (default: 'PBE1PBE/6-31G*').
    resp_theory : str, optional
        Theory for RESP calculation (default: 'HF/6-31G*').
    net_charge : float, optional
        Net charge for the molecule (default: 0.0).
    force_gaussian_rerun : bool, optional
        Whether to force rerun of Gaussian (default: False).
    minimize : bool, optional
        Whether to perform minimization (default: True).

    Attributes
    ----------
    in_mol2 : Path
        Path to the input mol2 file.
    out_gaussian_log : Path
        Path to the output Gaussian log file.
    opt_theory : str
        Theory for optimization.
    resp_theory : str
        Theory for RESP calculation.
    net_charge : float
        Net charge for the molecule.
    force_gaussian_rerun : bool
        Whether to force rerun of Gaussian.
    gaussian_cwd : Path
        Directory for Gaussian calculations.
    minimize : bool
        Whether to perform minimization.
    label : str
        Label for the calculation.
    """

    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_mol2 = Path(main_input)
        self.out_gaussian_log = Path(kwargs["out_gaussian_log"])

        self._validate_input_paths(**kwargs)
        self.opt_theory = kwargs.get("opt_theory", "PBE1PBE/6-31G*")
        self.resp_theory = kwargs.get("resp_theory", "HF/6-31G*")
        self.net_charge = kwargs.get("net_charge", 0.0)
        self.force_gaussian_rerun = kwargs.get("force_gaussian_rerun", False)
        self.gaussian_cwd = Path(self.cwd, "gaussianCalcs")
        self.minimize = kwargs.get("minimize", True)

        self.label = self.out_gaussian_log.stem

        return

    def _validate_input_paths(self, **kwargs):
        apply_gaussian_env_paths(self, kwargs)

    def setup(self, name_template: str) -> bool:
        """
        Set up Gaussian input and output files for the calculation.

        Parameters
        ----------
        name_template : str
            Template name for input/output files.

        Returns
        -------
        bool
            True if Gaussian calculation is already complete, False otherwise.
        """
        self.in_com = self.gaussian_cwd / f"{name_template}.com"
        self.out_log = self.gaussian_cwd / f"{name_template}.log"
        self._add_outputs(self.out_log)

        # __init__ tries to set up the coordinates object, but it may not have been available at init time.
        print(f"Setting up Gaussian calculations in {self.gaussian_cwd}")
        self.logger.info(f"Setting up Gaussian calculations in {self.gaussian_cwd}")
        if not getattr(self, "coord_object", None):
            self.coord_object = Coordinates(self.in_mol2, filetype="pdb")
        self.gaussian_cwd.mkdir(exist_ok=True)

        stageheader = _gaussian_link0_header(
            self.nproc, self.mem, chk=f"{self.in_mol2.stem}.antechamber.chk"
        )

        # Set up the Gaussian Block - it does not yet write anything,
        # so this part can be set up before the Gaussian calculations are run.
        gau = GaussianWriter(self.in_com)
        if self.minimize:
            n_atoms = len(self.coord_object.get_elements())
            opt_keyword = _gaussian_opt_keyword(n_atoms)
            self.logger.info(
                f"Gaussian optimization keyword: {opt_keyword} "
                f"({n_atoms} atoms; CalcFC if N < {_CALCFC_MAX_ATOMS})"
            )
            gau.add_block(
                GaussianInput(
                    command=f"#P {self.opt_theory} {opt_keyword}",
                    initial_coordinates=self.coord_object.get_coordinates(),
                    elements=self.coord_object.get_elements(),
                    charge=self.net_charge,
                    header=stageheader,
                )
            )
            gau.add_block(
                GaussianInput(
                    command=f"#P {self.resp_theory} GEOM(AllCheck) Guess(Read) NoSymm Pop=mk IOp(6/33=2) GFInput GFPrint",
                    charge=self.net_charge,
                    header=stageheader,
                )
            )
        else:
            extra = (
                "GEOM(AllCheck) Guess(Read) "
                if getattr(self, "_esp_from_chk", False)
                else ""
            )
            gau.add_block(
                GaussianInput(
                    command=(
                        f"#P {self.resp_theory} {extra}"
                        "NoSymm Pop=mk IOp(6/33=2) GFInput GFPrint"
                    ),
                    initial_coordinates=self.coord_object.get_coordinates(),
                    elements=self.coord_object.get_elements(),
                    charge=self.net_charge,
                    header=stageheader,
                )
            )

        gau_complete = _should_skip_gaussian_job(
            force_rerun=bool(self.force_gaussian_rerun),
            final_log=self.out_gaussian_log,
            cwd_log=self.out_log,
            logger=self.logger,
        )

        if not gau_complete:
            gau.write(dry_run=False)

        return gau_complete

    def _run(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        """Write the Gaussian input (unless skipped) and run the job."""
        if self.setup(self.label):
            self.logger.info(
                "Gaussian minimize/RESP already complete; skipping execution"
            )
            return
        _run_gaussian_com_and_promote(self, dry_run=dry_run)


class GaussianRESP(GaussianMinimizeRESP):
    """RESP-only Gaussian job that reads geometry and guess from the chk file.

    Same public constructor as historically: ``out_gaussian_log``,
    ``resp_theory``, ``net_charge``. Internally this is
    ``GaussianMinimizeRESP(minimize=False)`` with ``GEOM(AllCheck) Guess(Read)``.
    """

    _esp_from_chk = True

    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        kwargs.setdefault("minimize", False)
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)


class StageGaussianRotation(AbstractStage):
    """Rotate the ligand and run a Gaussian ESP job at each orientation.

    Supports two protocols:

    * ``so3_n28`` - 28 deterministic quaternion-packed SO(3) orientations
    * ``legacy_euler`` - historical Rx/Ry Euler grid (requires ``alpha``,
      ``beta``, ``gamma`` lists)

    Output logs are named ``{out_gaussian_label}_rot_*.log`` so
    :class:`~ligandparam.stages.Resp.StageMultiRespFit` can discover them
    regardless of protocol.

    Parameters
    ----------
    stage_name : str
        Stage name.
    main_input : path-like
        Input mol2 used for the rotated ESP jobs.
    cwd : path-like
        Working directory (Gaussian files go under ``cwd/gaussianCalcs``).
    out_gaussian_label : str
        Filename prefix for ``.com`` / ``.log`` outputs.
    orientation_protocol : {"legacy_euler", "so3_n28"}, optional
        Orientation generator. Default ``legacy_euler`` when used as a
        standalone stage; recipes such as FreeLigand typically pass ``so3_n28``.
    alpha, beta, gamma : list of float, optional
        Euler angles in degrees (Rx, Ry, Rz). Required for ``legacy_euler``.
    resp_theory : str, optional
        Theory for the ESP / RESP single-point jobs.
    net_charge : float, optional
        Net molecular charge.
    force_gaussian_rerun : bool, optional
        If False (default), skip orientation logs that already show
        ``Normal termination``. If True, rerun every orientation.
    nproc : int, optional
        Total core budget for this stage. Concurrent jobs and per-job
        ``%NProc`` are chosen so ``n_workers * %NProc <= nproc``.
    mem : int, optional
        Total memory budget (GB). Split across concurrent jobs so
        ``n_workers * %MEM <= mem`` (each job also gets at least 4 GB
        unless the allocation is smaller).
    """

    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_mol2 = Path(main_input)
        self.out_gaussian_label = kwargs["out_gaussian_label"]

        self._validate_input_paths(**kwargs)
        self.opt_theory = kwargs.get("opt_theory", "HF/6-31G*")
        self.resp_theory = kwargs.get("resp_theory", "HF/6-31G*")
        self.net_charge = kwargs.get("net_charge", 0.0)
        self.force_gaussian_rerun = kwargs.get("force_gaussian_rerun", False)
        self.gaussian_cwd = Path(self.cwd, "gaussianCalcs")

        self.orientation_protocol = kwargs.get("orientation_protocol", "legacy_euler")
        if self.orientation_protocol == "legacy_euler":
            if "alpha" not in kwargs or "beta" not in kwargs or "gamma" not in kwargs:
                raise ValueError(
                    "legacy_euler requires alpha, beta, and gamma angle lists"
                )
            self.alpha = [float(a) for a in kwargs["alpha"]]
            self.beta = [float(b) for b in kwargs["beta"]]
            self.gamma = [float(g) for g in kwargs["gamma"]]
        elif self.orientation_protocol == "so3_n28":
            self.alpha = []
            self.beta = []
            self.gamma = []
        else:
            raise ValueError(
                "orientation_protocol must be 'legacy_euler' or 'so3_n28'"
            )

        self.in_com_template = Path(self.gaussian_cwd, f"{self.out_gaussian_label}.com")
        self.xyz = Path(self.gaussian_cwd, f"{self.out_gaussian_label}_rotations.xyz")

    def _validate_input_paths(self, **kwargs):
        apply_gaussian_env_paths(self, kwargs)

    def _orientation_coordinates(self):
        """Yield stable filename suffixes and coordinates for each orientation."""
        if self.orientation_protocol == "so3_n28":
            quaternions = get_quaternion_pack("so3_n28")
            minimum_angle = minimum_pairwise_rotation_angle(quaternions)
            self.logger.info(
                f"Using {len(quaternions)}-point SO(3) quaternion pack "
                f"(minimum pairwise angle {minimum_angle:.2f} degrees)"
            )
            for index, quaternion in enumerate(quaternions):
                rotation = quaternion_to_matrix(quaternion)
                yield f"q{index:03d}", self.coord_object.rotate_matrix(rotation)
            return

        for alpha, beta, gamma in product(self.alpha, self.beta, self.gamma):
            suffix = f"{alpha:0.2f}_{beta:0.2f}_{gamma:0.2f}"
            yield suffix, self.coord_object.rotate(
                alpha=alpha, beta=beta, gamma=gamma
            )

    def _n_orientation_count(self) -> int:
        """Number of orientation jobs this stage will write."""
        if self.orientation_protocol == "so3_n28":
            return len(get_quaternion_pack("so3_n28"))
        return len(self.alpha) * len(self.beta) * len(self.gamma)

    def setup(self, name_template: str) -> bool:
        """
        Set up Gaussian input and output files for the rotation calculations.

        Parameters
        ----------
        name_template : str or Path
            Template name for input/output files. Accepted as either a string
            or a Path (callers may pass either).

        Returns
        -------
        bool
            Always returns False (rotation calculations are not pre-completed).
        """
        job_nproc = getattr(self, "_job_nproc", None) or self.nproc
        job_mem = getattr(self, "_job_mem", None) or self.mem
        self.header = _gaussian_link0_header(job_nproc, job_mem)

        # __init__ tries to set up the coordinates object, but it may not have been available at init time.
        if not getattr(self, "coord_object", None):
            self.coord_object = Coordinates(self.in_mol2, filetype="pdb")
        self.gaussian_cwd.mkdir(exist_ok=True)
        logger.info(f"Setting up Gaussian calculations in {self.gaussian_cwd}")
        print(f"Setting up Gaussian calculations in {self.gaussian_cwd}")

        # Some recipes pass a Path while others pass a string.
        name_label = Path(name_template).name

        store_coords = []
        self.in_coms = []
        self.out_logs = []
        elements = self.coord_object.get_elements()
        for orientation_suffix, test_rotation in self._orientation_coordinates():
            store_coords.append(test_rotation)
            # Keep "<rotation label>_*.log" stable: StageMultiRespFit discovers
            # these files by that prefix regardless of the orientation protocol.
            in_com = self.gaussian_cwd / f"{name_label}_rot_{orientation_suffix}.com"
            print(f"--> Writing Gaussian input file: {in_com}")
            self.in_coms.append(in_com)
            newgau = GaussianWriter(in_com)
            newgau.add_block(
                GaussianInput(
                    command=f"#P {self.resp_theory} SCF(Conver=6) NoSymm Test Pop=mk IOp(6/33=2) GFInput GFPrint",
                    initial_coordinates=test_rotation,
                    elements=elements,
                    charge=self.net_charge,
                    header=self.header,
                )
            )
            # Always write the Gaussian input file
            newgau.write(dry_run=False)

            out_log = self.gaussian_cwd / f"{name_label}_rot_{orientation_suffix}.log"
            self.out_logs.append(out_log)
            self._add_outputs(out_log)

        # Write the coordinates to a "trajectory" file
        self.write_rotation(store_coords, name_label)

        return False

    def _run(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        """Pool orientation ESP jobs. ``nproc`` / ``mem`` are node budgets."""
        import multiprocessing as mp

        from ligandparam.runtime.ProgressBoard import JobBoardWatcher, JobProgressStore

        n_orients = self._n_orientation_count()
        n_workers, job_nproc, job_mem = split_gaussian_orientation_budget(
            self.nproc, n_orients, self.mem
        )
        self._rotation_n_workers = n_workers
        self._job_nproc = job_nproc
        self._job_mem = job_mem
        self.logger.info(
            f"Gaussian rotation parallel plan: {n_orients} job(s), "
            f"nproc={self.nproc} mem={self.mem}GB -> "
            f"{n_workers} worker(s) x %NProc={job_nproc} x %MEM={job_mem}GB"
        )
        self.setup(self.out_gaussian_label)

        if self.force_gaussian_rerun:
            self.logger.info(
                "force_gaussian_rerun (-O): all orientation ESP jobs will be re-run"
            )

        status_path = self.gaussian_cwd / ".rot_progress.json"
        board_path = self.gaussian_cwd / "ROT_STATUS.txt"
        store = JobProgressStore(
            status_path,
            collection_key="orientations",
            id_header="Angle",
            title="Gaussian orientation ESP - live status",
            empty_hint="no orientations registered yet",
            detail_hint_label="Per-orientation Gaussian logs",
        )

        pending = []
        for in_com, out_log in zip(self.in_coms, self.out_logs):
            job_id = _orientation_id_from_paths(in_com, out_log)
            already_done = (
                not self.force_gaussian_rerun
                and _gaussian_log_is_complete(out_log)
            )
            if already_done:
                store.register(
                    job_id,
                    status="skipped",
                    stage="finished",
                    detail=f"already complete | {out_log.name}",
                    log_path=str(out_log),
                )
                self.logger.info(
                    "Skipping complete orientation ESP: %s", out_log.name
                )
                continue
            if out_log.exists() and not self.force_gaussian_rerun:
                self.logger.info(
                    "Incomplete orientation ESP log (%s); will re-run this job",
                    out_log.name,
                )
            store.register(
                job_id,
                status="queued",
                stage="queued",
                detail=f"{out_log.name} | %NProc={job_nproc} %MEM={job_mem}GB",
                log_path=str(out_log),
            )
            pending.append(
                {
                    "cwd": str(self.gaussian_cwd),
                    "in_com": in_com.name,
                    "out_log": out_log.name,
                    "job_id": job_id,
                    "status_path": str(status_path),
                    "job_nproc": int(job_nproc),
                    "job_mem": int(job_mem),
                    "force": bool(self.force_gaussian_rerun),
                    "dry_run": bool(dry_run),
                    "gaussian_root": self.gaussian_root,
                    "gauss_exedir": self.gauss_exedir,
                    "gaussian_binary": self.gaussian_binary,
                    "gaussian_scratch": self.gaussian_scratch,
                }
            )

        total = len(self.in_coms)
        finished = total - len(pending)
        self.logger.info(
            "Gaussian rotation status board: %s "
            "(%s already complete, %s pending)",
            board_path,
            finished,
            len(pending),
        )
        watcher = JobBoardWatcher(
            store,
            board_path=board_path,
            logger=self.logger,
            interval_sec=5.0,
            log_root_hint=str(self.gaussian_cwd / "*_rot_*.log"),
            thread_name="rot-progress-board",
        )
        watcher.start()
        try:
            if not pending:
                return

            workers = min(n_workers, len(pending))
            if dry_run or workers <= 1:
                for i, job in enumerate(pending):
                    _run_gaussian_rotation_job(job)
                    self._print_status(finished + i + 1, total)
                return

            ctx = mp.get_context("spawn")
            with ctx.Pool(processes=workers) as pool:
                for i, _result in enumerate(
                    pool.imap_unordered(_run_gaussian_rotation_job, pending)
                ):
                    self._print_status(finished + i + 1, total)
        finally:
            watcher.stop()

        return

    def _print_status(self, count, total_count):
        """Log progress through the orientation set."""
        percent = count / total_count * 100
        self.logger.info(f"Current Rotation Progress: {percent:.2f}%")

    def write_rotation(self, coords, name_template: str):
        """Write all rotated frames to ``{label}_rotations.xyz``."""
        self.logger.info(f"--> Writing rotations to file: gaussianCalcs/{name_template}_rotations.xyz")
        with open(self.xyz, "w") as file_obj:
            for frame in coords:
                SimpleXYZ(file_obj, frame)


class StageGaussianToMol2(AbstractStage):
    """
    Convert Gaussian output to mol2 format and assign charges to the mol2 file.

    Parameters
    ----------
    stage_name : str
        The name of the stage.
    main_input : Union[Path, str]
        Path to the input Gaussian log file.
    cwd : Union[Path, str]
        Current working directory.
    template_mol2 : str
        Path to the template mol2 file.
    out_mol2 : str
        Path to the output mol2 file.
    net_charge : float, optional
        Net charge for the molecule (default: 0.0).
    atom_type : str, optional
        Atom type (default: 'gaff2').
    force_gaussian_rerun : bool, optional
        Whether to force rerun of Gaussian (default: False).

    Attributes
    ----------
    in_log : Path
        Path to the input Gaussian log file.
    template_mol2 : Path
        Path to the template mol2 file.
    out_mol2 : Path
        Path to the output mol2 file.
    temp1_mol2 : Path
        Path to the first temporary mol2 file.
    temp2_mol2 : Path
        Path to the second temporary mol2 file.
    net_charge : float
        Net charge for the molecule.
    atom_type : str
        Atom type.
    force_gaussian_rerun : bool
        Whether to force rerun of Gaussian.
    gaussian_cwd : Path
        Directory for Gaussian calculations.
    """

    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_log = Path(main_input)
        self.template_mol2 = Path(kwargs["template_mol2"])
        self.out_mol2 = Path(kwargs["out_mol2"])
        self.temp1_mol2 = Path(self.cwd, f"{self.out_mol2.stem}.tmp1.mol2")
        self.temp2_mol2 = Path(self.cwd, f"{self.out_mol2.stem}.tmp2.mol2")
        self.net_charge = kwargs.get("net_charge", 0.0)
        self.atom_type = kwargs.get("atom_type", "gaff2")

        self._validate_input_paths(**kwargs)
        self.net_charge = kwargs.get("net_charge", 0.0)
        self.force_gaussian_rerun = kwargs.get("force_gaussian_rerun", False)
        self.gaussian_cwd = Path(self.cwd, "gaussianCalcs")

        self._add_outputs(self.out_mol2)
        self.add_required(self.in_log)

    def _validate_input_paths(self, **kwargs) -> None:
        apply_gaussian_env_paths(self, kwargs)

    def _run(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        """Convert the Gaussian log to mol2 and copy charges from the template."""
        warnings.filterwarnings("ignore")

        # Convert from gaussian to mol2
        ante = Antechamber(cwd=self.cwd, logger=self.logger, nproc=self.nproc)
        ante.call(i=self.in_log, fi="gout", o=self.temp1_mol2, fo="mol2", pf="y", at=self.atom_type, an="no", nc=self.net_charge, dry_run=dry_run)

        # Assign the charges
        if not dry_run:
            u1 = mda.Universe(self.template_mol2)
            u2 = mda.Universe(self.temp1_mol2)
            assert len(u1.atoms) == len(u2.atoms), "Number of atoms in the two files do not match"

            u2.atoms.charges = u1.atoms.charges
            Mol2Writer(u2, self.temp2_mol2, selection="all").write()

        # Use antechamber to clean up the mol2 format
        ante = Antechamber(cwd=self.cwd, logger=self.logger, nproc=self.nproc)
        ante.call(i=self.temp2_mol2, fi="mol2", o=self.out_mol2, fo="mol2", pf="y", at=self.atom_type, an="no", nc=self.net_charge, dry_run=dry_run)

        return


# Back-compat alias
StageGaussiantoMol2 = StageGaussianToMol2
