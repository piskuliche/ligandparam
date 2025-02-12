import logging
import os
from abc import abstractmethod
from typing_extensions import override
import subprocess
from pathlib import Path

from ligandparam.log import get_logger


class SimpleInterface:

    @abstractmethod
    def __init__(self, *args, **kwargs) -> None:
        """ This class is a simple interface to call external programs.

        This class is a simple interface to call external programs. It is designed to be subclassed
        and the method attribute set to the desired program. The call method will then call the program
        with the specified arguments.
        
        """
        pass

    def set_method(self, method):
        self.method = method
        return

    def call(self, **kwargs):
        dry_run = False
        if "dry_run" in kwargs:
            dry_run = kwargs['dry_run']
            del kwargs['dry_run']

        command = [self.method]
        shell = False
        for key, value in kwargs.items():
            if key == "inp_pipe":
                command.extend([f'<', str(value)])
                shell = True
            elif key == "out_pipe":
                command.extend([f'>', str(value)])
                shell = True
            else:
                if value is not None:
                    command.extend([f'-{key}', str(value)])

        if dry_run:
            self.logger.info(f"Command: {' '.join(command)}")
        else:
            env = os.environ
            if hasattr(self, "nproc"):
                # Prevent antechamber from using more threads than available
                env["OMP_NUM_THREADS"] = str(self.nproc)
            self.logger.info('\t' + ' '.join(command))
            p = subprocess.run(command, shell=shell, encoding='utf-8', cwd=self.cwd, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, env=env)
            if p.returncode != 0:
                self.logger.error(f"Command at {self.cwd} failed.")
                self.logger.error(p.stdout)
                self.logger.error(p.stderr)
                raise RuntimeError

        return


class Antechamber(SimpleInterface):
    @override
    def __init__(self, *args, **kwargs) -> None:
        """ This class is a simple interface to call the Antechamber program. """
        try:
            self.cwd = Path(kwargs["cwd"])
        except KeyError:
            raise ValueError(f"ERROR: missing `cwd` arg with a path to the workdir.")
        self.logger = kwargs.get('logger', get_logger())
        self.nproc = kwargs.get('nproc', 1)
        self.set_method('antechamber')
        return


class ParmChk(SimpleInterface):
    @override
    def __init__(self, *args, **kwargs) -> None:
        """ This class is a simple interface to call the ParmChk program. """
        try:
            self.cwd = Path(kwargs["cwd"])
        except KeyError:
            raise ValueError(f"ERROR: missing `cwd` arg with a path to the workdir.")
        self.logger = kwargs.get('logger', get_logger())
        self.set_method('parmchk2')
        return


class Leap(SimpleInterface):
    @override
    def __init__(self, *args, **kwargs) -> None:
        """ This class is a simple interface to call the Leap program. """
        try:
            self.cwd = Path(kwargs["cwd"])
        except KeyError:
            raise ValueError(f"ERROR: missing `cwd` arg with a path to the workdir.")
        self.logger = kwargs.get('logger', get_logger())
        self.set_method('tleap')
        return


class Gaussian(SimpleInterface):
    @override
    def __init__(self, *args, **kwargs) -> None:
        """ This class is a simple interface to call the Gaussian program. """
        try:
            self.cwd = Path(kwargs["cwd"])
        except KeyError:
            raise ValueError(f"ERROR: missing `cwd` arg with a path to the workdir.")
        for opt in ("gaussian_root", "gauss_exedir", "gaussian_binary", "gaussian_scratch"):
            try:
                setattr(self, opt, kwargs[opt])
            except KeyError:
                raise ValueError(f"ERROR: Please provide {opt} option as a keyword argument.")

        self.logger = kwargs.get('logger', get_logger())
        self.set_method(str(self.gaussian_binary))
        return

    def call(self, **kwargs):
        """ This function calls the Gaussian program with the specified arguments,
        however, it works slightly differently than the other interfaces. The Gaussian
        interface for some reason isn't compatible with the subprocess.run() function
        so we instead write a bash script to call the program and then execute the script."""

        dry_run = False
        if "dry_run" in kwargs:
            dry_run = kwargs['dry_run']
            del kwargs['dry_run']

        command = [self.method]
        shell = False
        for key, value in kwargs.items():
            if key == "inp_pipe":
                command.extend([f'<', str(value)])
                shell = True
            elif key == "out_pipe":
                command.extend([f'>', str(value)])
                shell = True
            else:
                if value is not None:
                    command.extend([f'-{key}', str(value)])

        self.write_bash(' '.join(command))
        bashcommand = 'bash temp_gaussian_sub.sh'

        if dry_run:
            self.logger.info(f"Command: {bashcommand}")
        else:
            self.logger.info('\t' + bashcommand)

            # Set the Gaussian environment variables if they weren't already set
            env = os.environ
            env.setdefault("g16root", str(self.gaussian_root))
            env.setdefault("GAUSS_EXEDIR", str(self.gauss_exedir))
            env.setdefault("GAUSS_SCRDIR", str(self.gaussian_scratch))

            p = subprocess.run(bashcommand, shell=shell, cwd=self.cwd, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, env=env)
            if p.returncode != 0:
                logging.error(f"Gaussian run at {self.cwd} failed.")
                logging.error(p.stdout)
                logging.error(p.stderr)
                raise RuntimeError

        return

    def write_bash(self, command):
        """ This function writes a bash script to call the Gaussian program
        with the specified arguments. """
        with open(self.cwd / "temp_gaussian_sub.sh", 'w') as f:
            f.write('#!/bin/bash\n\n')
            f.write(command)
            f.write('\n')
        return
