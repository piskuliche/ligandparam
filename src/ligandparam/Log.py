"""Logging for standalone ligandparam (same contract as GitLab ``main``).

This module is the only logging setup the package needs. It does not import
``runtime.Console`` or ALPS. ALPS can later wrap or parse these lines; the
formats below are the contract.

Stream (``logger="stream"`` / :func:`set_stream_logger`)
    The record message only, on stdout. No timestamp, no level, no tag.
    Prints the ligandparam ASCII banner + a random quote once.

File (``logger="file"`` / :func:`set_file_logger`)
    ``YYYY-mm-dd HH:MM:SS - LEVEL - message``
    Banner goes to stdout (unprefixed); the quote is also written to the file.

CLI file logs (``lig-getparam``) may insert the package version after the
level, matching the historical ``ligandparam_getparam`` formatter::

    YYYY-mm-dd HH:MM:SS - LEVEL - 1.6.1 message
"""

from __future__ import annotations

import logging
import os
import random
import sys
from pathlib import Path
from typing import Optional

from . import __logging_name__, __version__

_QUOTES_FILE = Path(__file__).resolve().parent / "pkgdata" / "quotes.txt"
_BANNER_PRINTED = False

# Historical ligandparam ASCII wordmark (standalone package, not ALPS).
_LOGO = r"""
.____    .__                         ._____________
|    |   |__| _________    ____    __| _/\______   \_____ ____________    _____
|    |   |  |/ ___\__  \  /    \  / __ |  |     ___/\__  \\_  __ \__  \  /     \
|    |___|  / /_/  > __ \|   |  \/ /_/ |  |    |     / __ \|  | \// __ \|  Y Y  \
|_______ \__\___  (____  /___|  /\____ |  |____|    (____  /__|  (____  /__|_|  /
        \/ /_____/     \/     \/      \/                 \/           \/      \/
""".strip(
    "\n"
)

_AUTHORS = (
    "Zeke Piskulich <piskulichz@gmail.com>",
    "German P. Barletta",
    "Timothy J. Giese",
    "Nate Levinzon <ndlevinzon@gmail.com>",
)

FILE_DATEFMT = "%Y-%m-%d %H:%M:%S"
FILE_FORMAT = "{asctime} - {levelname} - {message}"
FILE_FORMAT_WITH_VERSION = "{asctime} - {levelname} - {version} {message}"


def default_quotes_path() -> Path:
    """Packaged quotes file (one quote per line)."""
    return _QUOTES_FILE


def load_quotes(path: Path | None = None) -> list[str]:
    """Read non-empty, non-comment lines from the quotes file."""
    quotes_path = Path(path) if path is not None else default_quotes_path()
    if not quotes_path.is_file():
        return []
    out: list[str] = []
    try:
        text = quotes_path.read_text(encoding="utf-8")
    except OSError:
        return []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if len(line) >= 2 and line[0] == line[-1] and line[0] in {'"', "'"}:
            line = line[1:-1].strip()
        if line:
            out.append(line)
    return out


def pick_quote(quotes_path: Path | None = None) -> Optional[str]:
    """Return one random quote, or None if the quotes file is empty."""
    quotes = load_quotes(quotes_path)
    if not quotes:
        return None
    return random.choice(quotes)


def format_reminder_line(quote: str) -> str:
    """Return the success reminder line for one quote."""
    return f"LIGANDPARAM reminds you: {quote}"


def format_startup_banner(
    *,
    version: str | None = None,
    quote: str | None = None,
) -> str:
    """ASCII wordmark + version + authors + optional quote (no log prefix)."""
    ver = version if version is not None else __version__
    authors = "\n".join(f"  {a}" for a in _AUTHORS)
    lines = [
        _LOGO,
        "",
        f"  ligandparam  v{ver}",
        "  Amber ligand parameterization",
        "",
        "  Authors:",
        authors,
    ]
    if quote:
        lines.extend(["", f"  {format_reminder_line(quote)}"])
    return "\n".join(lines) + "\n"


def print_startup_banner(
    *,
    logger: logging.Logger | None = None,
    stream=None,
    force: bool = False,
    quotes_path: Path | None = None,
) -> bool:
    """Print the ligandparam banner once (stdout, unprefixed).

    A random quote is included in the banner. When ``logger`` is given, the
    same reminder line is also written through the logger so file logs keep
    it. Spawn workers skip a reprint via ``LIGANDPARAM_BANNER_PRINTED``.

    Returns
    -------
    bool
        ``True`` if the banner was printed on this call.
    """
    global _BANNER_PRINTED
    if not force:
        if _BANNER_PRINTED or os.environ.get("LIGANDPARAM_BANNER_PRINTED"):
            return False
    quote = pick_quote(quotes_path)
    text = format_startup_banner(quote=quote)
    out = stream if stream is not None else sys.stdout
    try:
        out.write(text if text.endswith("\n") else text + "\n")
        out.flush()
    except OSError:
        return False
    if logger is not None and quote:
        logger.info(format_reminder_line(quote))
    _BANNER_PRINTED = True
    os.environ["LIGANDPARAM_BANNER_PRINTED"] = "1"
    return True


def log_success_quote(
    logger: logging.Logger | None = None,
    *,
    quotes_path: Path | None = None,
) -> Optional[str]:
    """Pick a random quote and emit it through ``logger`` (or stdout).

    Uses :meth:`logging.Logger.info` when ``logger`` is given so the line
    follows the same file/stream format as the rest of the run. No-op when
    the quotes file is missing or empty.
    """
    quote = pick_quote(quotes_path)
    if not quote:
        return None
    line = format_reminder_line(quote)
    if logger is not None:
        logger.info(line)
    else:
        print(line, flush=True)
    return quote


def dihed_correct_ok(result, *, dry_run: bool = False) -> bool:
    """True when a dihedral-correct result wrote a frcmod with no failed jobs."""
    if dry_run or not isinstance(result, dict):
        return False
    out = result.get("merged_frcmod") or result.get("out_frcmod")
    if not out or not Path(out).is_file():
        return False
    for frag in result.get("fragments") or []:
        if not isinstance(frag, dict):
            continue
        status = str(frag.get("status") or "").lower()
        if status == "failed" or frag.get("error"):
            return False
    report = result.get("merge_report") or {}
    if isinstance(report, dict) and report.get("errors"):
        return False
    return True


def _drop_null_handlers(logger: logging.Logger) -> None:
    logger.handlers = [
        h for h in logger.handlers if not isinstance(h, logging.NullHandler)
    ]


def file_formatter(*, include_version: bool = False) -> logging.Formatter:
    """Return the GitLab-main file formatter (optional CLI version field)."""
    if include_version:
        return logging.Formatter(
            FILE_FORMAT_WITH_VERSION,
            style="{",
            datefmt=FILE_DATEFMT,
            defaults={"version": __version__},
        )
    return logging.Formatter(
        FILE_FORMAT,
        style="{",
        datefmt=FILE_DATEFMT,
    )


def get_logger() -> logging.Logger:
    """Return the package logger with a null handler until configured."""
    logger = logging.getLogger(__logging_name__)
    logger.setLevel(logging.INFO)
    if not logger.handlers:
        logger.addHandler(logging.NullHandler())
    return logger


def set_stream_logger(logging_level: int = logging.INFO) -> logging.Logger:
    """Log messages to stdout with no extra prefix (GitLab ``main`` behavior)."""
    logger = logging.getLogger(__logging_name__)
    logger.setLevel(logging_level)
    _drop_null_handlers(logger)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging_level)
    logger.addHandler(stream_handler)
    print_startup_banner()
    return logger


def set_file_logger(
    logfilename: Path,
    logname: str = None,
    filemode: str = "a",
    *,
    include_version: bool = False,
) -> logging.Logger:
    """Log to a file using the GitLab ``main`` ``asctime - LEVEL - message`` format.

    Prints the ASCII banner to stdout. The quote is also written to the file.
    """
    if logname is None:
        logname = __logging_name__
    logger = logging.getLogger(logname)
    logger.setLevel(logging.INFO)
    _drop_null_handlers(logger)
    file_handler = logging.FileHandler(filename=logfilename, mode=filemode)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(file_formatter(include_version=include_version))
    logger.addHandler(file_handler)
    print_startup_banner(logger=logger)
    return logger
