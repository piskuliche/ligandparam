"""
This module provides classes for parametrizing ligands and managing recipes.

Classes
-------
Parametrization
    A class for parametrizing ligands using various stages.

Recipe
    A subclass of Parametrization for managing ligand recipes.
"""

import logging
from pathlib import Path
from typing import Optional, Union
from typing_extensions import override

from ligandparam.driver import Driver
from ligandparam.log import get_logger, set_stream_logger, set_file_logger


class Parametrization(Driver):
    """
    A class for parametrizing ligands using various stages.

    Parameters
    ----------
    in_filename : Union[Path, str]
        The input filename of the ligand.
    cwd : Union[Path, str]
        The current working directory.
    *args : tuple
        Additional positional arguments.
    **kwargs : dict
        Additional keyword arguments.

    Keyword Args
    ------------
    label : str, optional
        A label for the ligand, by default the stem of `in_filename`.
    leaprc : list, optional
        A list of leaprc files to use, by default ["leaprc.gaff2"].
    logger : Union[str, logging.Logger], optional
        The logger to use. Can be "file", "stream", or a logging.Logger instance.

    Attributes
    ----------
    in_filename : Path
        The resolved path to the input file.
    label : str
        The label for the ligand.
    cwd : Path
        The current working directory.
    stages : list
        A list of stages to run.
    leaprc : list
        A list of leaprc files to use.
    logger : logging.Logger
        The logger instance.

    Raises
    ------
    ValueError
        If an invalid logger type is provided.
    """

    @override
    def __init__(self, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs):
        """
        The rough approach to using this class is to generate a new Parametrization class, and then generate self.stages as a list
        of stages that you want to run.

        Parameters
        ----------
        in_filename : Union[Path, str]
            The input filename of the ligand.
        cwd : Union[Path, str]
            The current working directory.
        *args : tuple
            Additional positional arguments.
        **kwargs : dict
            Additional keyword arguments.

        Keyword Args
        ------------
        label : str, optional
            A label for the ligand, by default the stem of `in_filename`.
        leaprc : list, optional
            A list of leaprc files to use, by default ["leaprc.gaff2"].
        logger : Union[str, logging.Logger], optional
            The logger to use. Can be "file", "stream", or a logging.Logger instance.

        Raises
        ------
        ValueError
            If an invalid logger type is provided.
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
        """
        Add a leaprc file to the list of leaprc files.

        Parameters
        ----------
        leaprc : str
            The name of the leaprc file to add.
        """
        self.leaprc.append(leaprc)


class Recipe(Parametrization):
    """
    A subclass of Parametrization for managing ligand recipes.
    """
    pass
