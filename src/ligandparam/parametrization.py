import logging
from pathlib import Path
from typing import Optional,  Union
from typing_extensions import override

from ligandparam.driver import Driver
from ligandparam.log import get_logger, set_stream_logger, set_file_logger


class Parametrization(Driver):
    @override
    def __init__(self, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs):
        """
        The rough approach to using this class is to generate a new Parametrization class, and then generate self.stages as a list
        of stages that you want to run.
        Args:
            in_filename (str): The in_filename of the ligand.
            cwd (Union[Path, str]): The current working directory.
            *args: Additional positional arguments.
            **kwargs: Additional keyword arguments.
        Keyword Args:
            name (str): The base name for the ligand.
            inputoptions (dict): A dictionary of input options, which should include 'name' or 'pdb_filename'.
        Raises:
            ValueError: If neither 'name' nor 'pdb_filename' is provided in inputoptions.
        """
        self.in_filename = Path(in_filename).resolve()
        self.label = kwargs.get("label", self.in_filename.stem)
        self.cwd = Path(cwd)
        self.stages = []
        self.leaprc = kwargs.get("leaprc", ["leaprc.gaff2"])
        try:
            logger = kwargs.pop("logger")
            if isinstance(logger, str):
                if logger == "file":
                    self.logger = set_file_logger(self.cwd / f"{self.label}.log")
                elif logger == "stream":
                    self.logger = set_stream_logger()
                else:
                    raise ValueError("Invalid input string for logger. Must be either 'file' or 'stream'.")
            elif isinstance(logger, logging.Logger):
                self.logger = logger
            else:
                raise ValueError("logger must be a string or a logging.Logger instance.")
        except KeyError:
            self.logger = get_logger()

    def add_leaprc(self, leaprc) -> None:
        self.leaprc.append(leaprc)


class Recipe(Parametrization):
    pass
