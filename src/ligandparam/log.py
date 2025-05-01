import sys
import logging
from pathlib import Path

from . import __logging_name__


def get_logger() -> logging.Logger:
    logger = logging.getLogger(__logging_name__)
    logger.setLevel(logging.INFO)
    logger.addHandler(logging.NullHandler())
    return logger

def set_stream_logger(logging_level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger(__logging_name__)
    logger.setLevel(logging_level)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging_level)
    logger.addHandler(stream_handler)

    return logger

def set_file_logger(logfilename: Path, logname: str = None, filemode: str = 'a') -> logging.Logger:
    if logname is None:
        # logname = Path(logfilename).stem
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
