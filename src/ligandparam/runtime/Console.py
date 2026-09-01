"""Console formatting helpers for Slurm-friendly stdout/stderr logging.

Format (single timestamp, hierarchical brackets)::

    YYYY-mm-dd HH:MM:SS [ligandparam] INFO: message
    YYYY-mm-dd HH:MM:SS [ligandparam] [gaussian] INFO: message

Leading ``[scope]`` tokens in the message body are peeled into the bracket
hierarchy so callers can keep writing ``log.info("[frag-twist] ...")`` without
duplicating tags. Console handlers always write to the real process streams
(``sys.__stdout__`` / ``sys.__stderr__``) so teeing stdout for fragment logs
does not double-prefix already-formatted logger lines.

All console writes are forced to ASCII (``+/-``, ``deg``, ``chi^2``, ...) so
latin-1 Slurm ``.out`` files do not mojibake.
"""

from __future__ import annotations

import errno
import logging
import re
import sys
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Sequence, TextIO

# Lines that already carry our console prefix (logger -> TeeTextIO).
_PREFIXED_LINE = re.compile(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \[")
# Leading ``[scope]`` tokens embedded in the message / print body.
_LEADING_SCOPE = re.compile(r"^\[([^\]]+)\]\s*")

_BANNER_PRINTED = False

# Readable ASCII stand-ins for symbols that otherwise mojibake on latin-1
# Slurm ``.out`` files. Anything else is ``backslashreplace``'d.
_ASCII_REPLACEMENTS = (
    ("\u00b1", "+/-"),  # plus-minus
    ("\u2213", "-/+"),
    ("\u00b0", " deg"),  # degree
    ("\u03c7", "chi"),
    ("\u03a7", "chi"),
    ("\u00b2", "^2"),
    ("\u00b3", "^3"),
    ("\u0394", "d"),
    ("\u03b4", "d"),
    ("\u2192", "->"),
    ("\u2190", "<-"),
    ("\u2194", "<->"),
    ("\u2014", "-"),
    ("\u2013", "-"),
    ("\u2212", "-"),
    ("\u00d7", "x"),
    ("\u2264", "<="),
    ("\u2265", ">="),
    ("\u2260", "!="),
    ("\u2248", "~="),
    ("\u221d", "propto"),
    ("\u00b7", "*"),
    ("\u2022", "*"),
    ("\u03bc", "u"),
    ("\u00b5", "u"),
    ("\u00c5", "Ang"),
    ("\u2026", "..."),
    ("\u2018", "'"),
    ("\u2019", "'"),
    ("\u201c", '"'),
    ("\u201d", '"'),
    ("\u00a0", " "),
)


def is_stale_file_handle(exc: BaseException) -> bool:
    """True for NFS/VAST ``ESTALE`` (Linux errno 116)."""
    if not isinstance(exc, OSError):
        return False
    estale = getattr(errno, "ESTALE", 116)
    if getattr(exc, "errno", None) in (estale, 116):
        return True
    msg = str(exc).lower()
    return "stale file handle" in msg or "estale" in msg


def _reopen_filehandler_stream(handler) -> bool:
    """Reopen a ``FileHandler`` stream after ESTALE. Return True on success."""
    opener = getattr(handler, "_open", None)
    if not callable(opener):
        return False
    try:
        stream = getattr(handler, "stream", None)
        if stream is not None:
            try:
                stream.close()
            except OSError:
                pass
        handler.stream = opener()
        return True
    except OSError:
        return False


def install_stale_handle_logging_guard() -> None:
    """Do not dump Python ``Logging error`` tracebacks on VAST stale fds.

    geomeTRIC's ``RawFileHandler`` writes a log on scratch. When that fd goes
    ESTALE, ``StreamHandler.flush`` raises and ``Handler.handleError`` prints
    a full traceback per optimizer step, flooding Slurm ``.out``. Swallow the
    error and reopen file-backed handlers so the opt can finish.
    """
    if getattr(logging.Handler.handleError, "_ffpopt_estale", False):
        return

    _orig_handleError = logging.Handler.handleError
    _orig_flush = logging.StreamHandler.flush

    def _handleError(self, record):
        exc = sys.exc_info()[1]
        if is_stale_file_handle(exc):
            _reopen_filehandler_stream(self)
            return
        return _orig_handleError(self, record)

    def _flush(self):
        try:
            _orig_flush(self)
        except OSError as exc:
            if not is_stale_file_handle(exc):
                raise
            _reopen_filehandler_stream(self)

    _handleError._ffpopt_estale = True  # type: ignore[attr-defined]
    logging.Handler.handleError = _handleError
    logging.StreamHandler.flush = _flush


def install_ase_futurewarning_filter() -> None:
    """Silence ASE ``ignore_bad_restart_file`` FutureWarning (stderr flood).

    ``ase.calculators.amber.SANDER`` still calls ``Calculator`` with extra
    positional args / that deprecated keyword. Each spawn worker would
    otherwise print the same warning on every calculator build.

    ``PYTHONWARNINGS`` is also set so multiprocessing spawn workers and
    geomeTRIC child interpreters apply the same filter at startup.
    """
    import os
    import warnings

    if getattr(install_ase_futurewarning_filter, "_done", False):
        return
    warnings.filterwarnings(
        "ignore",
        category=FutureWarning,
        message=r".*ignore_bad_restart_file.*",
        module=r"ase\.calculators(\.calculator)?",
    )
    warnings.filterwarnings(
        "ignore",
        category=FutureWarning,
        message=r".*ignore_bad_restart_file.*",
    )
    extra = "ignore:.*ignore_bad_restart_file:FutureWarning"
    existing = os.environ.get("PYTHONWARNINGS", "")
    if "ignore_bad_restart_file" not in existing:
        os.environ["PYTHONWARNINGS"] = (
            f"{existing},{extra}" if existing else extra
        )
    install_ase_futurewarning_filter._done = True  # type: ignore[attr-defined]


def ascii_for_stdio(text: str) -> str:
    """Return ``text`` encoded as ASCII suitable for stdout / Slurm logs."""
    if not text or text.isascii():
        return text
    out = text
    for src, dst in _ASCII_REPLACEMENTS:
        if src in out:
            out = out.replace(src, dst)
    if out.isascii():
        return out
    return out.encode("ascii", errors="backslashreplace").decode("ascii")


class _AsciiStdio:
    """Proxy that forces ``write`` / ``writelines`` through :func:`ascii_for_stdio`."""

    def __init__(self, inner: TextIO) -> None:
        self._inner = inner
        self._lp_ascii_stdio = True

    def write(self, data: str) -> int:
        if not isinstance(data, str):
            data = str(data)
        try:
            return self._inner.write(ascii_for_stdio(data))
        except OSError as exc:
            if is_stale_file_handle(exc):
                return 0
            raise

    def flush(self):
        try:
            return self._inner.flush()
        except OSError as exc:
            if is_stale_file_handle(exc):
                return None
            raise

    def writelines(self, lines) -> None:
        for line in lines:
            self.write(line)

    def __getattr__(self, name: str):
        return getattr(self._inner, name)


def ensure_ascii_stdio() -> None:
    """Wrap ``sys.stdout`` / ``sys.stderr`` so print() cannot emit non-ASCII."""
    install_stale_handle_logging_guard()
    install_ase_futurewarning_filter()
    for attr in ("stdout", "stderr"):
        stream = getattr(sys, attr, None)
        if stream is None or getattr(stream, "_lp_ascii_stdio", False):
            continue
        setattr(sys, attr, _AsciiStdio(stream))


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


def package_version() -> str:
    """Return the installed ``ligandparam`` version string."""
    try:
        from ligandparam import __version__

        return str(__version__)
    except Exception:
        try:
            from importlib.metadata import version

            return version("ligandparam")
        except Exception:
            return "unknown"


def format_startup_banner(*, version: str | None = None) -> str:
    """ASCII logo + version + authors (delegates to :mod:`ligandparam.Log`)."""
    from ligandparam.Log import format_startup_banner as _format

    return _format(version=version)


def format_run_banner(
    title: str,
    fields: Sequence[tuple[str, str]],
    *,
    width: int = 72,
) -> str:
    """Fixed-width ASCII card for the current job (no timestamp prefix).

    Used under the startup logo so Slurm ``.out`` files open with a greppable
    summary (mode, ligand, nproc, log paths) before wavefront spam.
    """
    width = max(52, min(int(width), 96))
    inner = width - 2
    title = ascii_for_stdio(str(title or "")).strip() or "RUN"
    if len(title) > inner - 2:
        title = title[: inner - 2]

    def _bar() -> str:
        return "+" + "-" * inner + "+"

    def _fit(line: str) -> str:
        if not line:
            return _bar()
        cap = line[0] if line[0] in "+|" else "|"
        end = line[-1] if line[-1] in "+|" else "|"
        body = line[1:-1] if len(line) >= 2 else ""
        if len(body) < inner:
            body = body + " " * (inner - len(body))
        elif len(body) > inner:
            body = body[:inner]
        return cap + body + end

    label_w = 8
    if fields:
        label_w = max(len(str(k)) for k, _ in fields)
        label_w = min(max(label_w, 8), 12)
    value_w = inner - 4 - label_w
    if value_w < 16:
        value_w = 16
        label_w = max(6, inner - 4 - value_w)

    lines = [_bar(), "|" + title.center(inner) + "|", _bar()]
    for key, raw in fields:
        label = ascii_for_stdio(str(key)).strip()
        val = ascii_for_stdio(str(raw if raw is not None else "")).replace("\n", " ").strip()
        chunks: list[str] = []
        text = val
        while text:
            chunks.append(text[:value_w])
            text = text[value_w:]
        if not chunks:
            chunks = [""]
        for i, chunk in enumerate(chunks):
            lab = label if i == 0 else ""
            lines.append(f"|  {lab:<{label_w}} {chunk:<{value_w}} |")
    lines.append(_bar())
    return "\n".join(_fit(line) for line in lines) + "\n"


def print_run_banner(text: str, *, stream: TextIO | None = None) -> None:
    """Write a run card to the real stdout (not through a fragment/batch tee)."""
    ensure_ascii_stdio()
    out = stream if stream is not None else sys.__stdout__
    body = ascii_for_stdio(text)
    if not body.endswith("\n"):
        body += "\n"
    try:
        out.write(body if body.startswith("\n") else "\n" + body)
        out.flush()
    except OSError:
        pass


def format_whole_ligand_run_banner(
    *,
    ligand: str,
    model: str = "qdpi2",
    nproc: int = 1,
    delta: int = 10,
    n_bonds: int = 0,
    n_batches: int | None = None,
    extras: str = "",
    work_dir: str = "",
) -> str:
    """ASCII card for ``--whole-ligand`` / AFFDO-style twist runs."""
    if n_batches is None:
        bonds_line = f"{int(n_bonds)} rotatable"
    else:
        bonds_line = (
            f"{int(n_bonds)} rotatable  ->  {int(n_batches)} sequential "
            f"batch(es)"
        )
    extra = extras.strip() if extras else "(none)"
    return format_run_banner(
        "WHOLE-LIGAND TWIST",
        [
            ("ligand", ligand),
            ("model", f"{model}    nproc={int(nproc)}    delta={int(delta)} deg"),
            ("bonds", bonds_line),
            ("extras", extra),
            ("work", work_dir),
            ("status", "WHOLE_STATUS.txt"),
            ("logs", "torsion_batch_XX/whole-twist.log"),
        ],
    )


def format_fragmented_run_banner(
    *,
    ligand: str,
    model: str = "qdpi2",
    nproc: int = 1,
    n_fragments: int = 0,
    work_dir: str = "",
) -> str:
    """ASCII card for the default scission / fragment twist path."""
    return format_run_banner(
        "FRAGMENTED TWIST",
        [
            ("ligand", ligand),
            ("model", f"{model}    nproc={int(nproc)}"),
            ("frags", f"{int(n_fragments)} fragment(s)"),
            ("work", work_dir),
            ("status", "FRAG_STATUS.txt"),
            ("logs", "<fragment>/frag-twist.log"),
        ],
    )


def print_startup_banner(
    *,
    stream: TextIO | None = None,
    force: bool = False,
) -> bool:
    """Print the project banner once for this job (stdout top).

    Uses a process-local flag **and** ``LIGANDPARAM_BANNER_PRINTED`` so spawned
    fragment / wavefront workers that re-import this module do not reprint.
    Do not call from ``attach_console_handlers`` - only from top-level CLI
    entry points.

    Returns
    -------
    bool
        ``True`` if the banner was printed on this call.
    """
    import os

    global _BANNER_PRINTED
    if not force:
        if _BANNER_PRINTED or os.environ.get("LIGANDPARAM_BANNER_PRINTED"):
            return False
    ensure_ascii_stdio()
    out = stream if stream is not None else sys.__stdout__
    from ligandparam.Log import print_startup_banner as _print

    ok = bool(_print(stream=out, force=True))
    if not ok:
        return False
    _BANNER_PRINTED = True
    os.environ["LIGANDPARAM_BANNER_PRINTED"] = "1"
    return True


def console_timestamp() -> str:
    """Return a local wall-clock timestamp suitable for console lines."""
    return time.strftime("%Y-%m-%d %H:%M:%S")


def _normalize_tags(tag: str | Sequence[str] | None = None, *extra: str) -> list[str]:
    tags: list[str] = []
    if tag is None:
        pass
    elif isinstance(tag, str):
        if tag:
            tags.append(tag)
    else:
        tags.extend(str(t) for t in tag if t)
    tags.extend(str(t) for t in extra if t)
    # Preserve order, drop duplicates.
    seen: set[str] = set()
    out: list[str] = []
    for t in tags:
        if t not in seen:
            seen.add(t)
            out.append(t)
    return out


def peel_leading_scopes(message: str) -> tuple[list[str], str]:
    """Split leading ``[scope]`` tokens from ``message``."""
    scopes: list[str] = []
    rest = message
    while True:
        match = _LEADING_SCOPE.match(rest)
        if not match:
            break
        scopes.append(match.group(1))
        rest = rest[match.end() :]
    return scopes, rest


def format_bracket_tags(*tags: str) -> str:
    """Return ``[a] [b]`` (no trailing space); empty if no tags."""
    return " ".join(f"[{t}]" for t in tags if t)


def format_console_line(
    message: str,
    *,
    tag: str | Sequence[str] | None = None,
    tags: Sequence[str] | None = None,
) -> str:
    """Prefix ``message`` once: ``TIMESTAMP [tag...] [peeled...] rest``."""
    if _PREFIXED_LINE.match(message.lstrip("\r")):
        # Already formatted (e.g. logging handler wrote into a TeeTextIO).
        lined = message if message.endswith("\n") else f"{message}\n"
        return ascii_for_stdio(lined)

    base = _normalize_tags(tag)
    if tags:
        base = _normalize_tags(base + list(tags))

    stamp = console_timestamp()
    text = message if message.endswith("\n") else f"{message}\n"
    if text == "\n":
        bracket = format_bracket_tags(*base)
        prefix = f"{stamp} {bracket} " if bracket else f"{stamp} "
        return prefix + "\n"

    parts = text.splitlines(keepends=True)
    out: list[str] = []
    for part in parts:
        newline = part.endswith("\n")
        body = part[:-1] if newline else part
        scopes, rest = peel_leading_scopes(body)
        bracket = format_bracket_tags(*_normalize_tags(base + scopes))
        prefix = f"{stamp} {bracket} " if bracket else f"{stamp} "
        line = prefix + rest
        out.append(line + ("\n" if newline else "\n"))
    return ascii_for_stdio("".join(out))


class TeeTextIO:
    """Write plain text to a log file and mirror tagged lines to a console stream."""

    def __init__(
        self,
        file_stream: TextIO,
        console_stream: TextIO,
        *,
        tag: str | Sequence[str],
        log_path: Path | None = None,
    ) -> None:
        self.file_stream = file_stream
        self.console_stream = console_stream
        self.tags = _normalize_tags(tag)
        self.log_path = Path(log_path) if log_path is not None else None
        self.encoding = getattr(file_stream, "encoding", "utf-8") or "utf-8"
        self._buf = ""

    def _reopen_log(self) -> bool:
        if self.log_path is None:
            return False
        try:
            try:
                self.file_stream.close()
            except OSError:
                pass
            self.file_stream = open(
                self.log_path, "a", encoding="utf-8", buffering=1
            )
            return True
        except OSError:
            return False

    def write(self, data: str) -> int:
        if not data:
            return 0
        data = ascii_for_stdio(data)
        try:
            self.file_stream.write(data)
        except OSError as exc:
            if not is_stale_file_handle(exc):
                raise
            if self._reopen_log():
                try:
                    self.file_stream.write(data)
                except OSError:
                    pass
        self._buf += data
        while "\n" in self._buf:
            line, self._buf = self._buf.split("\n", 1)
            try:
                self.console_stream.write(
                    format_console_line(line + "\n", tags=self.tags)
                )
            except OSError as exc:
                if not is_stale_file_handle(exc):
                    raise
        return len(data)

    def flush(self) -> None:
        if self._buf:
            try:
                self.console_stream.write(
                    format_console_line(self._buf, tags=self.tags)
                )
            except OSError as exc:
                if not is_stale_file_handle(exc):
                    raise
            self._buf = ""
        try:
            self.file_stream.flush()
        except OSError as exc:
            if is_stale_file_handle(exc):
                self._reopen_log()
        try:
            self.console_stream.flush()
        except OSError:
            pass

    def isatty(self) -> bool:
        return False

    def fileno(self) -> int:
        return self.file_stream.fileno()

    def writable(self) -> bool:
        return True

    def close(self) -> None:
        self.flush()


class _MaxLevelFilter(logging.Filter):
    """Allow records with ``levelno <= max_level``."""

    def __init__(self, max_level: int) -> None:
        super().__init__()
        self.max_level = max_level

    def filter(self, record: logging.LogRecord) -> bool:
        return record.levelno <= self.max_level


class HierarchicalConsoleFormatter(logging.Formatter):
    """``TIMESTAMP [base...] [peeled...] LEVEL: message`` with one timestamp."""

    def __init__(self, *tags: str, datefmt: str = "%Y-%m-%d %H:%M:%S") -> None:
        super().__init__(datefmt=datefmt)
        self.base_tags = _normalize_tags(tags)

    def format(self, record: logging.LogRecord) -> str:
        scopes, rest = peel_leading_scopes(record.getMessage())
        bracket = format_bracket_tags(*_normalize_tags(self.base_tags + scopes))
        asctime = self.formatTime(record, self.datefmt)
        level = record.levelname
        if bracket:
            head = f"{asctime} {bracket} {level}: {rest}"
        else:
            head = f"{asctime} {level}: {rest}"
        if record.exc_info:
            if not record.exc_text:
                record.exc_text = self.formatException(record.exc_info)
        if record.exc_text:
            if head[-1:] != "\n":
                head += "\n"
            head += record.exc_text
        return ascii_for_stdio(head)


def console_formatter(
    tag: str | Sequence[str] | None = None,
    *extra_tags: str,
) -> logging.Formatter:
    """Formatter: ``YYYY-mm-dd HH:MM:SS [tag...] LEVEL: message``."""
    tags = _normalize_tags(tag, *extra_tags)
    return HierarchicalConsoleFormatter(*tags)


def attach_console_handlers(
    logger: logging.Logger,
    *,
    tag: str | Sequence[str] = "ligandparam",
    level: int = logging.INFO,
) -> None:
    """Attach stdout (INFO) and stderr (WARNING+) handlers if not already present.

    Writes to ``sys.__stdout__`` / ``sys.__stderr__`` so temporary stdio tees
    cannot double-prefix logger output.
    """
    for handler in logger.handlers:
        if getattr(handler, "_lp_console_marker", None):
            return

    ensure_ascii_stdio()
    tags = _normalize_tags(tag)
    marker = "ligandparam.console:" + "/".join(tags)
    fmt = console_formatter(tags)

    out = logging.StreamHandler(_AsciiStdio(sys.__stdout__))
    out.setLevel(level)
    out.addFilter(_MaxLevelFilter(logging.INFO))
    out.setFormatter(fmt)
    out._lp_console_marker = marker  # type: ignore[attr-defined]
    logger.addHandler(out)

    err = logging.StreamHandler(_AsciiStdio(sys.__stderr__))
    err.setLevel(logging.WARNING)
    err.setFormatter(fmt)
    err._lp_console_marker = marker  # type: ignore[attr-defined]
    logger.addHandler(err)

    logger.setLevel(min(logger.level or level, level))
    # Avoid duplicate lines if a parent also has console handlers.
    logger.propagate = False


@contextmanager
def tee_stdio_to_file(
    log_path: Path,
    *,
    tag: str | Sequence[str] = "ffpopt",
) -> Iterator[None]:
    """Write stdout/stderr to ``log_path`` and mirror tagged lines to the console.

    Useful for wavefront ``print`` spam under parallel fragment pools: the
    per-fragment ``.log`` stays complete, and Slurm ``.out`` / ``.err`` still
    receive the same content with a single timestamp and hierarchical tags.
    """
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = open(log_path, "a", encoding="utf-8", buffering=1)
    old_out, old_err = sys.stdout, sys.stderr
    # Mirror to the streams we are replacing (real console, or a test capture).
    # Logger console handlers use sys.__stdout__/__stderr__, so they are not
    # double-prefixed by this tee.
    tee_out = TeeTextIO(fh, old_out, tag=tag, log_path=log_path)
    tee_err = TeeTextIO(fh, old_err, tag=tag, log_path=log_path)
    try:
        sys.stdout = tee_out  # type: ignore[assignment]
        sys.stderr = tee_err  # type: ignore[assignment]
        yield
    finally:
        try:
            tee_out.flush()
            tee_err.flush()
        except OSError:
            pass
        sys.stdout = old_out
        sys.stderr = old_err
        try:
            fh.close()
        except OSError:
            pass
