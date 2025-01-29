import logging
from . import __logging_name__


def get_logger() -> logging.Logger:
    # If they ever change this loggerDict, we're screwed
    if __logging_name__ in logging.Logger.manager.loggerDict.keys():
        logger = logging.getLogger(__logging_name__)
    else:
        logger = logging.getLogger(__logging_name__)
        logger.setLevel(logging.INFO)
        logger.addHandler(logging.NullHandler())
    return logger