import logging
from . import __logging_name__


def get_logger() -> logging.Logger:
    logger = logging.getLogger(__logging_name__)
    logger.setLevel(logging.INFO)
    logger.addHandler(logging.NullHandler())
    return logger
