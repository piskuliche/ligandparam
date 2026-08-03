"""Tests for handler management in ligandparam.log."""

import logging

import pytest

from ligandparam import __logging_name__
from ligandparam.log import get_logger, set_file_logger, set_stream_logger


@pytest.fixture(autouse=True)
def clean_logger():
    """Each test starts from a logger with no handlers and restores them afterwards."""
    logger = logging.getLogger(__logging_name__)
    saved = list(logger.handlers)
    logger.handlers.clear()
    yield logger
    for handler in logger.handlers:
        handler.close()
    logger.handlers.clear()
    logger.handlers.extend(saved)


def test_get_logger_does_not_accumulate_handlers(clean_logger):
    """get_logger() is the default for every stage and interface, so repeated calls
    used to leave one NullHandler per stage on the shared logger."""
    for _ in range(10):
        get_logger()
    assert len(clean_logger.handlers) == 1


def test_set_stream_logger_does_not_duplicate(clean_logger):
    """A second recipe in the same process used to print every line twice."""
    set_stream_logger()
    set_stream_logger()
    assert len(clean_logger.handlers) == 1


def test_set_file_logger_does_not_duplicate(tmp_path, clean_logger):
    log = tmp_path / "run.log"
    set_file_logger(log)
    set_file_logger(log)
    assert len(clean_logger.handlers) == 1


def test_distinct_log_files_each_get_a_handler(tmp_path, clean_logger):
    set_file_logger(tmp_path / "a.log")
    set_file_logger(tmp_path / "b.log")
    assert len(clean_logger.handlers) == 2


def test_messages_are_written_once(tmp_path, clean_logger):
    log = tmp_path / "run.log"
    logger = set_file_logger(log)
    set_file_logger(log)
    logger.info("hello")
    for handler in clean_logger.handlers:
        handler.flush()
    assert log.read_text().count("hello") == 1
