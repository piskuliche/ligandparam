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
                self.logger.error(f"Command at {self.cwd} failed.")
                self.logger.error(p.stdout)
                self.logger.error(p.stderr)
                raise RuntimeError(p.stderr)

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
        """This function calls the Gaussian program with the specified arguments,
        however, it works slightly differently than the other interfaces. The Gaussian
        interface for some reason isn't compatible with the subprocess.run() function
        so we instead write a bash script to call the program and then execute the script."""

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

        self.write_bash(" ".join(command))
        bashcommand = "bash temp_gaussian_sub.sh"

        if dry_run:
            self.logger.info(f"Command: {bashcommand}")
        else:
            self.logger.info("\t" + bashcommand)

            # Set the Gaussian environment variables if they weren't already set
            env = self.set_environment()

            p = subprocess.run(
                bashcommand, shell=shell, cwd=self.cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env
            )
            if p.returncode != 0:
                self.logger.error(f"Gaussian run at {self.cwd} failed.")
                self.logger.error(p.stdout)
                self.logger.error(p.stderr)
                raise RuntimeError

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
        env = os.environ
        if not env.get("g16root") and self.gaussian_root:
            env["g16root"] = str(self.gaussian_root)
        if not env.get("GAUSS_EXEDIR") and self.gauss_exedir:
            env["GAUSS_EXEDIR"] = str(self.gauss_exedir)
        if not env.get("GAUSS_SCRDIR") and self.gaussian_scratch:
            env["GAUSS_SCRDIR"] = str(self.gaussian_scratch)
        return env
