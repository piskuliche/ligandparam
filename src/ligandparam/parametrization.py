from pathlib import Path
from typing import Union
from typing_extensions import override

from ligandparam.driver import Driver
from ligandparam.log import get_logger


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
        self.in_filename = Path(in_filename)
        self.label = kwargs.get("label", self.in_filename.stem)
        self.cwd = Path(cwd)
        self.logger = kwargs.get("logger", get_logger())
        self.stages = []
        self.leaprc = kwargs.get("leaprc", ["leaprc.gaff2"])

    def add_leaprc(self, leaprc) -> None:
        self.leaprc.append(leaprc)


class Recipe(Parametrization):
    pass
