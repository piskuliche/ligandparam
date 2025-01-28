from pathlib import Path
from typing import Union

from ligandparam.driver import Driver
from ligandparam.io.coordinates import Coordinates
from ligandparam.stages import *


class Parametrization(Driver):
    def __init__(self, name: str, cwd: Union[Path, str], *args, **kwargs):
        """
        The rough approach to using this class is to generate a new Parametrization class, and then generate self.stages as a list
        of stages that you want to run.
        Args:
            name (str): The name of the ligand.
            cwd (Union[Path, str]): The current working directory.
            *args: Additional positional arguments.
            **kwargs: Additional keyword arguments.
        Keyword Args:
            name (str): The base name for the ligand.
            inputoptions (dict): A dictionary of input options, which should include 'name' or 'pdb_filename'.
        Raises:
            ValueError: If neither 'name' nor 'pdb_filename' is provided in inputoptions.
        """

        self.name = name
        self.cwd = Path(cwd)
        self.stages = []
        self.leaprc = []
        if hasattr(kwargs, 'leaprc'):
            self.leaprc = kwargs['leaprc']

    def add_leaprc(self, leaprc) -> None:
        self.leaprc.append(leaprc)


class Recipe(Parametrization):
    def __init__(self, name: str, cwd: Union[Path, str], *args, **kwargs):
        super().__init__(name, cwd, *args, **kwargs)
        return
