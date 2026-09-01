"""Wrappers that invoke external Amber/Gaussian tools via subprocess."""

import logging
import os
from abc import abstractmethod
from typing_extensions import override
import subprocess
from pathlib import Path

from ligandparam.Log import get_logger


def _pipe_text(blob) -> str:
    """Decode subprocess PIPE bytes for logs / exception messages."""
    if blob is None:
        return ""
    if isinstance(blob, bytes):
        return blob.decode("utf-8", errors="replace")
    return str(blob)


def _subprocess_failure_message(cwd, p, tool: str = "Command") -> str:
    """Build a RuntimeError body from a failed ``subprocess.run`` result."""
    stdout = _pipe_text(p.stdout).strip()
    stderr = _pipe_text(p.stderr).strip()
    rc = p.returncode
    hint = ""
    if rc in (-9, 9, 137):
        hint = (
            " Process was SIGKILL'd (often the Slurm OOM killer). "
            "For orientation ESP, --mem is the node budget and is split "
            "across concurrent jobs."
        )

    def _tail(text: str, n: int = 1500) -> str:
        return text[-n:] if len(text) > n else text

    parts = [f"{tool} at {cwd} failed (returncode={rc}).{hint}"]
    if stderr:
        parts.append(f"stderr: {_tail(stderr)}")
    if stdout:
        parts.append(f"stdout: {_tail(stdout)}")
    return " ".join(parts)


class SimpleInterface:
    """Base wrapper for calling an external program.

    Subclasses set ``method`` to the executable name. :meth:`call` runs the
    program in ``cwd`` with optional logging and dry-run support.

    Attributes
    ----------
    method : str
        Executable name or path.
    logger : logging.Logger
        Logger used for command output.
    cwd : Path
        Working directory for the subprocess.
    nproc : int
        Processor hint forwarded by callers that support it.
    """

    @abstractmethod
    def __init__(self, *args, **kwargs) -> None:
        """
        Initialize the SimpleInterface class.

        This class is designed to be subclassed, with the `method` attribute set to the desired program. The `call` method will then execute the program with the specified arguments.
        """
        pass

    def set_method(self, method):
        """
        Set the method to call the external program.

        Parameters
        ----------
        method : str
            The name of the external program to call.
        """
        self.method = method
        return

    def call(self, **kwargs):
        """
        Call the external program with the specified arguments.

        Parameters
        ----------
        **kwargs : dict
            Keyword arguments to pass to the external program. Special keys include:
            - `dry_run` (bool): If True, log the command without executing it.
            - `inp_pipe` (str): Input file to pipe into the program.
            - `out_pipe` (str): Output file to pipe the program's output.

        Raises
        ------
        RuntimeError
            If the external program returns a non-zero exit code.
        """
        dry_run = False
        if "dry_run" in kwargs:
            dry_run = kwargs["dry_run"]
            del kwargs["dry_run"]

        command = [self.method]
        shell = False
        for key, value in kwargs.items():
            if key == "inp_pipe":
                command.extend(["<", str(value)])
                shell = True
            elif key == "out_pipe":
                command.extend([">", str(value)])
                shell = True
            else:
                if value is not None:
                    command.extend([f"-{key}", str(value)])

        if dry_run:
            self.logger.info(f"Command: {' '.join(command)}")
        else:
            env = os.environ
            if hasattr(self, "nproc"):
                # Prevent antechamber from using more threads than available
                env["OMP_NUM_THREADS"] = str(self.nproc)
            self.logger.info("\t" + " ".join(command))
            p = subprocess.run(
                command,
                shell=shell,
                encoding="utf-8",
                cwd=self.cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=env,
            )
            if p.returncode != 0:
                msg = _subprocess_failure_message(self.cwd, p, tool=self.method)
                self.logger.error(msg)
                raise RuntimeError(msg)

        return


class Antechamber(SimpleInterface):
    """
    Interface to call the Antechamber program.

    This class provides a simple interface to execute the Antechamber program.

    Parameters
    ----------
    *args : list
        Additional arguments to pass to the interface.
    **kwargs : dict
        Additional keyword arguments. Must include:
        - `cwd` (str): Path to the working directory.

    Attributes
    ----------
    cwd : Path
        The current working directory to run the program in.
    logger : logging.Logger
        The logger to use for logging.
    nproc : int
        The number of processors to use for the program.
    """

    @override
    def __init__(self, *args, **kwargs) -> None:
        """
        Initialize the Antechamber interface.

        Parameters
        ----------
        *args : list
            Additional arguments to pass to the interface.
        **kwargs : dict
            Additional keyword arguments. Must include:
            - `cwd` (str): Path to the working directory.

        Raises
        ------
        ValueError
            If the `cwd` argument is missing.
        """
        try:
            self.cwd = Path(kwargs["cwd"])
        except KeyError:
            raise ValueError(f"ERROR: missing `cwd` arg with a path to the workdir.")

        self.logger = kwargs.get("logger", get_logger())
        self.nproc = kwargs.get("nproc", 1)
        self.set_method("antechamber")
        return


class ParmChk(SimpleInterface):
    @override
    def __init__(self, *args, **kwargs) -> None:
        """This class is a simple interface to call the ParmChk program."""
        try:
            self.cwd = Path(kwargs["cwd"])
        except KeyError:
            raise ValueError(f"ERROR: missing `cwd` arg with a path to the workdir.")
        self.logger = kwargs.get("logger", get_logger())
        self.set_method("parmchk2")
        return


class Leap(SimpleInterface):
    @override
    def __init__(self, *args, **kwargs) -> None:
        """This class is a simple interface to call the Leap program."""
        try:
            self.cwd = Path(kwargs["cwd"])
        except KeyError:
            raise ValueError(f"ERROR: missing `cwd` arg with a path to the workdir.")
        self.logger = kwargs.get("logger", get_logger())
        self.set_method("tleap")
        return


class Gaussian(SimpleInterface):
    @override
    def __init__(self, *args, **kwargs) -> None:
        """This class is a simple interface to call the Gaussian program."""
        try:
            self.cwd = Path(kwargs["cwd"])
        except KeyError:
            raise ValueError(f"ERROR: missing `cwd` arg with a path to the workdir.")
        for opt in ("gaussian_root", "gauss_exedir", "gaussian_binary", "gaussian_scratch"):
            try:
                setattr(self, opt, kwargs.get(opt, ""))
            except KeyError:
                raise ValueError(f"ERROR: Please provide {opt} option as a keyword argument.")

        self.logger = kwargs.get("logger", get_logger())
        self.set_method(str(self.gaussian_binary))
        return

    def call(self, **kwargs):
        """Call Gaussian via a per-job bash wrapper (subprocess-safe).

        Extra kwargs (stripped before building the Gaussian command):

        - ``dry_run`` (bool): log the bash command without running it
        - ``script_name`` (str): bash wrapper filename under ``cwd`` (unique
          per concurrent job; default ``_gau_<stem>.sh`` from ``inp_pipe``)
        - ``scratch`` (path-like): ``GAUSS_SCRDIR`` for this job (default
          ``cwd/tmp/scratch_<stem>`` so concurrent jobs do not collide)
        """
        dry_run = bool(kwargs.pop("dry_run", False))
        script_name = kwargs.pop("script_name", None)
        scratch = kwargs.pop("scratch", None)

        command = [self.method]
        shell = False
        for key, value in kwargs.items():
            if key == "inp_pipe":
                command.extend(["<", str(value)])
                shell = True
            elif key == "out_pipe":
                command.extend([">", str(value)])
                shell = True
            else:
                if value is not None:
                    command.extend([f"-{key}", str(value)])

        inp = kwargs.get("inp_pipe")
        stem = Path(str(inp)).stem if inp else "job"
        if script_name is None:
            script_name = f"_gau_{stem}.sh"
        if scratch is None:
            if inp:
                scratch = self.cwd / "tmp" / f"scratch_{stem}"
            elif self.gaussian_scratch:
                scratch = self.gaussian_scratch

        self.write_bash(" ".join(command), script_name=script_name)
        bashcommand = f"bash {script_name}"

        if dry_run:
            self.logger.info(f"Command: {bashcommand}")
        else:
            self.logger.info("\t" + bashcommand)
            env = self.set_environment(scratch=scratch)
            p = subprocess.run(
                bashcommand,
                shell=shell,
                cwd=self.cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=env,
            )
            if p.returncode != 0:
                msg = _subprocess_failure_message(self.cwd, p, tool="Gaussian")
                self.logger.error(msg)
                raise RuntimeError(msg)

        return

    def write_bash(self, command, script_name="temp_gaussian_sub.sh"):
        """Write a bash script that invokes Gaussian with the given command line."""
        path = self.cwd / script_name
        with open(path, "w") as f:
            f.write("#!/bin/bash\n\n")
            f.write(command)
            f.write("\n")
        return path

    def set_environment(self, scratch=None) -> dict:
        """Build a subprocess env with Gaussian paths; never mutate ``os.environ``.

        Parameters
        ----------
        scratch : path-like, optional
            Job-specific ``GAUSS_SCRDIR``. Falls back to ``gaussian_scratch``,
            then any existing ``GAUSS_SCRDIR`` in the parent environment.
        """
        env = os.environ.copy()
        if self.gaussian_root:
            env["g16root"] = str(self.gaussian_root)
        if self.gauss_exedir:
            env["GAUSS_EXEDIR"] = str(self.gauss_exedir)
        scratch_dir = scratch if scratch is not None else self.gaussian_scratch
        if scratch_dir:
            Path(scratch_dir).mkdir(parents=True, exist_ok=True)
            env["GAUSS_SCRDIR"] = str(scratch_dir)
        return env
