"""
This module provides logging utilities for the ligandparam package.

Functions
---------
get_logger() -> logging.Logger
    Returns a logger with a null handler.

set_stream_logger(logging_level: int = logging.INFO) -> logging.Logger
    Sets up a logger to output to the standard output stream.

set_file_logger(logfilename: Path, logname: str = None, filemode: str = 'a') -> logging.Logger
    Sets up a logger to output to a file.
"""

import sys
import logging
from pathlib import Path

from . import __logging_name__


def get_logger() -> logging.Logger:
    """
    Get a logger with a null handler.

    Returns
    -------
    logging.Logger
        A logger instance with a null handler.
    """
    logger = logging.getLogger(__logging_name__)
    logger.setLevel(logging.INFO)
    logger.addHandler(logging.NullHandler())
    return logger

def set_stream_logger(logging_level: int = logging.INFO) -> logging.Logger:
    """
    Set up a logger to output to the standard output stream.

    Parameters
    ----------
    logging_level : int, optional
        The logging level to set for the logger, by default logging.INFO.

    Returns
    -------
    logging.Logger
        A logger instance configured to output to the standard output stream.
    """
    logger = logging.getLogger(__logging_name__)
    logger.setLevel(logging_level)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging_level)
    logger.addHandler(stream_handler)

    return logger

def set_file_logger(logfilename: Path, logname: str = None, filemode: str = 'a') -> logging.Logger:
    """
    Set up a logger to output to a file.

    Parameters
    ----------
    logfilename : Path
        The path to the log file.
    logname : str, optional
        The name of the logger, by default None. If None, the module's logging name is used.
    filemode : str, optional
        The mode to open the log file, by default 'a'.

    Returns
    -------
    logging.Logger
        A logger instance configured to output to the specified file.
    """
    if logname is None:
        logname = __logging_name__
    logger = logging.getLogger(logname)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    file_handler = logging.FileHandler(filename=logfilename, mode=filemode)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger
