import contextlib
import logging
import os
from abc import abstractmethod
from typing_extensions import override
import subprocess
from pathlib import Path

from ligandparam.log import get_logger


class SimpleInterface:
    """
    A simple interface to call external programs.

    This class is designed to be subclassed, with the `method` attribute set to the desired program. The `call` method will then execute the program with the specified arguments.

    Parameters
    ----------
    *args : list
        Additional arguments to pass to the subclass.
    **kwargs : dict
        Additional keyword arguments to pass to the subclass.

    Attributes
    ----------
    method : str
        The method to call the external program.
    logger : logging.Logger
        The logger to use for logging.
    cwd : Path
        The current working directory to run the program in.
    nproc : int
        The number of processors to use for the program.
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

        # Redirections are handled with real file handles rather than by appending
        # "<"/">" to the argument list. Passing a list with shell=True runs only
        # argv[0] and silently discards every argument, which meant piped calls
        # executed the bare binary and then reported success.
        command = [self.method]
        inp_pipe = None
        out_pipe = None
        for key, value in kwargs.items():
            if key == "inp_pipe":
                inp_pipe = value
            elif key == "out_pipe":
                out_pipe = value
            else:
                if value is not None:
                    command.extend([f"-{key}", str(value)])

        display = " ".join(command)
        if inp_pipe is not None:
            display += f" < {inp_pipe}"
        if out_pipe is not None:
            display += f" > {out_pipe}"

        if dry_run:
            self.logger.info(f"Command: {display}")
            return

        env = os.environ.copy()
        if hasattr(self, "nproc"):
            # Prevent antechamber from using more threads than available
            env["OMP_NUM_THREADS"] = str(self.nproc)
        self.logger.info("\t" + display)

        with contextlib.ExitStack() as stack:
            stdin = stdout = None
            if inp_pipe is not None:
                stdin = stack.enter_context(open(self._resolve(inp_pipe), "r"))
            if out_pipe is not None:
                stdout = stack.enter_context(open(self._resolve(out_pipe), "w"))
            p = subprocess.run(
                command,
                encoding="utf-8",
                cwd=self.cwd,
                stdin=stdin,
                stdout=stdout if stdout is not None else subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=env,
            )
        if p.returncode != 0:
            self.logger.error(f"Command at {self.cwd} failed.")
            if p.stdout:
                self.logger.error(p.stdout)
            if p.stderr:
                self.logger.error(p.stderr)
            raise RuntimeError(p.stderr)

        return

    def _resolve(self, filename) -> Path:
        """Resolve a redirection target relative to the interface's working directory.

        The shell used to interpret "<"/">" relative to ``cwd``; keep that behaviour
        now that the redirections are opened directly.
        """
        path = Path(filename)
        return path if path.is_absolute() else Path(self.cwd, path)


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
            setattr(self, opt, kwargs.get(opt))

        if not self.gaussian_binary:
            self.gaussian_binary = "g16"

        self.logger = kwargs.get("logger", get_logger())
        self.set_method(str(self.gaussian_binary))
        return

    def call(self, **kwargs):
        """This function calls the Gaussian program with the specified arguments,
        however, it works slightly differently than the other interfaces. The Gaussian
        interface for some reason isn't compatible with the subprocess.run() function
        so we instead write a bash script to call the program and then execute the script."""

        dry_run = False
        if "dry_run" in kwargs:
            dry_run = kwargs["dry_run"]
            del kwargs["dry_run"]

        if not self.method:
            raise ValueError(
                "ERROR: no Gaussian binary configured. Pass `gaussian_binary` (e.g. 'g16') "
                "as a keyword argument.")

        command = [self.method]
        for key, value in kwargs.items():
            # The redirections stay as shell text here: this command is written into a
            # bash script rather than executed directly.
            if key == "inp_pipe":
                command.extend(["<", str(value)])
            elif key == "out_pipe":
                command.extend([">", str(value)])
            else:
                if value is not None:
                    command.extend([f"-{key}", str(value)])

        self.write_bash(" ".join(command))
        bashcommand = ["bash", "temp_gaussian_sub.sh"]

        if dry_run:
            self.logger.info(f"Command: {' '.join(bashcommand)}")
        else:
            self.logger.info("\t" + " ".join(bashcommand))

            # Set the Gaussian environment variables if they weren't already set
            env = self.set_environment()

            p = subprocess.run(
                bashcommand, cwd=self.cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env
            )
            if p.returncode != 0:
                self.logger.error(f"Gaussian run at {self.cwd} failed.")
                self.logger.error(p.stdout)
                self.logger.error(p.stderr)
                raise RuntimeError(f"Gaussian run at {self.cwd} failed.")

        return

    def write_bash(self, command):
        """This function writes a bash script to call the Gaussian program
        with the specified arguments."""
        with open(self.cwd / "temp_gaussian_sub.sh", "w") as f:
            f.write("#!/bin/bash\n\n")
            f.write(command)
            f.write("\n")
        return

    def set_environment(self) -> dict:
        # Copy: assigning into os.environ directly would leak these settings into the
        # parent process and every later stage.
        env = os.environ.copy()
        if not env.get("g16root") and self.gaussian_root:
            env["g16root"] = str(self.gaussian_root)
        if not env.get("GAUSS_EXEDIR") and self.gauss_exedir:
            env["GAUSS_EXEDIR"] = str(self.gauss_exedir)
        if not env.get("GAUSS_SCRDIR") and self.gaussian_scratch:
            env["GAUSS_SCRDIR"] = str(self.gaussian_scratch)
        return env
