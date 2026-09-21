"""Shared ASCII job-status board for parallel workers.

Workers write stage updates into a JSON status file. A parent process (or
watcher thread) renders a fixed-width table so interleaved job output does not
obscure which unit of work is active.
"""

from __future__ import annotations

import json
import os
import re
import time
from pathlib import Path
from typing import Any, Mapping, Optional, TextIO

PathLike = str | Path


def atomic_write_text(path: Path, text: str) -> None:
    """Write ``text`` to ``path`` via a same-directory temp file + ``os.replace``."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


# Back-compat private alias
_atomic_write_text = atomic_write_text


class DirLock:
    """Exclusive lock via ``mkdir`` (works on NFS / Windows without extras)."""

    def __init__(self, lock_dir: Path, *, timeout_sec: float = 30.0) -> None:
        self.lock_dir = Path(lock_dir)
        self.timeout_sec = timeout_sec

    def __enter__(self) -> None:
        deadline = time.monotonic() + self.timeout_sec
        while True:
            try:
                os.mkdir(self.lock_dir)
                return
            except FileExistsError:
                if time.monotonic() >= deadline:
                    try:
                        os.rmdir(self.lock_dir)
                    except OSError:
                        pass
                    time.sleep(0.05)
                    continue

    def __exit__(self, *exc: object) -> None:
        try:
            os.rmdir(self.lock_dir)
        except OSError:
            pass


# Back-compat private alias
_DirLock = DirLock


class JobProgressStore:
    """Process-shared job status table backed by a JSON file."""

    def __init__(
        self,
        path: PathLike,
        *,
        collection_key: str = "jobs",
        id_header: str = "Job",
        title: str = "Live job status",
        empty_hint: str = "no jobs registered yet",
        detail_hint_label: str = "Detail logs",
    ) -> None:
        self.path = Path(path)
        self.lock_dir = self.path.with_suffix(self.path.suffix + ".lock")
        self.collection_key = collection_key
        self.id_header = id_header
        self.title = title
        self.empty_hint = empty_hint
        self.detail_hint_label = detail_hint_label

    def _read_unlocked(self) -> dict[str, Any]:
        if not self.path.is_file():
            return {self.collection_key: {}}
        try:
            data = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return {self.collection_key: {}}
        if not isinstance(data, dict):
            return {self.collection_key: {}}
        items = data.get(self.collection_key)
        if not isinstance(items, dict):
            data[self.collection_key] = {}
        return data

    def snapshot(self) -> dict[str, dict[str, Any]]:
        """Return ``{job_id: status_dict}``."""
        with _DirLock(self.lock_dir):
            data = self._read_unlocked()
        return dict(data.get(self.collection_key) or {})

    def register(
        self,
        job_id: str,
        *,
        status: str = "queued",
        stage: str = "queued",
        detail: str = "",
        log_path: str | None = None,
        **extra: Any,
    ) -> None:
        """Mark a job as present before workers start."""
        self.update(
            job_id,
            status=status,
            stage=stage,
            detail=detail,
            log_path=log_path,
            **extra,
        )

    def update(
        self,
        job_id: str,
        *,
        status: Optional[str] = None,
        stage: Optional[str] = None,
        detail: Optional[str] = None,
        error: Optional[str] = None,
        log_path: Optional[str] = None,
        **extra: Any,
    ) -> None:
        """Merge fields for one job and rewrite the JSON store."""
        with _DirLock(self.lock_dir):
            data = self._read_unlocked()
            items = data.setdefault(self.collection_key, {})
            entry = dict(items.get(job_id) or {})
            prev_status = str(entry.get("status") or "")
            prev_stage = str(entry.get("stage") or "")
            entry["id"] = job_id
            now = time.time()
            entry["updated"] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(now))
            if status is not None:
                entry["status"] = str(status)
            if stage is not None:
                entry["stage"] = str(stage)
            new_status = str(entry.get("status") or "")
            new_stage = str(entry.get("stage") or "")
            if new_status == "running":
                status_changed = status is not None and prev_status != "running"
                stage_changed = stage is not None and new_stage != prev_stage
                if status_changed or stage_changed or not entry.get("started_epoch"):
                    entry["started"] = entry["updated"]
                    entry["started_epoch"] = now
            elif new_status == "queued":
                entry.pop("started", None)
                entry.pop("started_epoch", None)
            if detail is not None:
                entry["detail"] = str(detail)
            if error is not None:
                entry["error"] = str(error)
            if log_path is not None:
                entry["log_path"] = str(log_path)
            for key, value in extra.items():
                if value is not None:
                    entry[key] = value
            items[job_id] = entry
            _atomic_write_text(self.path, json.dumps(data, indent=2) + "\n")

    def render_board(self, *, log_root_hint: str | None = None) -> str:
        """Return an ASCII status table for the current snapshot."""
        jobs = self.snapshot()
        return format_job_board(
            jobs,
            title=self.title,
            id_header=self.id_header,
            empty_hint=self.empty_hint,
            detail_hint_label=self.detail_hint_label,
            log_root_hint=log_root_hint,
        )


def job_id_sort_key(job_id: str) -> tuple:
    """Natural sort so ``fragment_2`` precedes ``fragment_10``."""
    parts: list[tuple[int, int | str]] = []
    for chunk in re.split(r"(\d+)", str(job_id)):
        if not chunk:
            continue
        if chunk.isdigit():
            parts.append((1, int(chunk)))
        else:
            parts.append((0, chunk.casefold()))
    return tuple(parts)


def jobs_fingerprint(jobs: Mapping[str, Mapping[str, Any]]) -> str:
    """Stable status signature (ignores elapsed / refresh timestamps)."""
    parts: list[str] = []
    for jid in sorted(jobs, key=job_id_sort_key):
        entry = jobs[jid] or {}
        parts.append(
            f"{jid}\t{entry.get('status', '')}\t{entry.get('stage', '')}\t"
            f"{entry.get('error', '')}"
        )
    return "\n".join(parts)


def elapsed_phrase(
    started: str | None,
    started_epoch: float | None = None,
) -> str:
    """Human elapsed time from a start stamp (prefer Unix epoch seconds)."""
    secs: int | None = None
    if started_epoch is not None:
        try:
            secs = max(0, int(time.time() - float(started_epoch)))
        except (TypeError, ValueError):
            secs = None
    if secs is None:
        if not started:
            return ""
        try:
            t0 = time.mktime(time.strptime(str(started), "%Y-%m-%d %H:%M:%S"))
        except (ValueError, OverflowError, OSError):
            return ""
        secs = max(0, int(time.time() - t0))
    if secs < 60:
        return f"{secs}s"
    mins, rem = divmod(secs, 60)
    if mins < 60:
        return f"{mins}m {rem}s"
    hours, rem_m = divmod(mins, 60)
    return f"{hours}h {rem_m}m"


def format_job_board(
    jobs: Mapping[str, Mapping[str, Any]],
    *,
    title: str = "Live job status",
    id_header: str = "Job",
    empty_hint: str = "no jobs registered yet",
    detail_hint_label: str = "Detail logs",
    log_root_hint: str | None = None,
) -> str:
    """Format ``{id: {status, stage, detail, ...}}`` as a fixed-width board."""
    ids = sorted(jobs.keys(), key=job_id_sort_key)
    col_id = max([len(id_header)] + [len(i) for i in ids] + [10])
    col_st = max(
        [len("Status")] + [len(str(jobs[i].get("status", ""))) for i in ids] + [8]
    )
    col_sg = max(
        [len("Stage")] + [len(str(jobs[i].get("stage", ""))) for i in ids] + [8]
    )
    col_dt = 52

    def row(jid: str, status: str, stage: str, detail: str) -> str:
        d = (detail or "")[: col_dt - 1]
        return (
            f" {jid:<{col_id}}  {status:<{col_st}}  {stage:<{col_sg}}  {d:<{col_dt}}"
        )

    width = col_id + col_st + col_sg + col_dt + 8
    bar = "=" * width
    sep = "-" * width
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    lines = [
        bar,
        f" {title}",
        f" refreshed {stamp}",
        bar,
        row(id_header, "Status", "Stage", "Detail"),
        sep,
    ]

    counts = {
        "done": 0,
        "running": 0,
        "queued": 0,
        "skipped": 0,
        "failed": 0,
        "other": 0,
    }
    for jid in ids:
        e = jobs[jid]
        status = str(e.get("status") or "?")
        stage = str(e.get("stage") or "-")
        detail = str(e.get("detail") or "")
        if status == "running":
            elapsed = elapsed_phrase(e.get("started"), e.get("started_epoch"))
            if elapsed:
                detail = (detail + " | " if detail else "") + f"elapsed {elapsed}"
        if e.get("error") and status == "failed":
            detail = (detail + " | " if detail else "") + str(e["error"])[:40]
        lines.append(row(jid, status, stage, detail))
        key = status if status in counts else "other"
        counts[key] = counts.get(key, 0) + 1

    if not ids:
        lines.append(row("(none)", "-", "-", empty_hint))

    lines.append(sep)
    summary = (
        f" {counts['done']} done | {counts['running']} running | "
        f"{counts['queued']} queued | {counts['skipped']} skipped | "
        f"{counts['failed']} failed"
    )
    if counts["other"]:
        summary += f" | {counts['other']} other"
    lines.append(summary)
    if log_root_hint:
        lines.append(f" {detail_hint_label}: {log_root_hint}")
    lines.append(bar)
    return "\n".join(lines) + "\n"


class JobBoardWatcher:
    """Background thread that refreshes a status board file and stdout.

    The JSON/board files update every ``interval_sec``. Stdout (or the logger
    fallback) reprints when a job's status fingerprint changes, on a
    heartbeat while something is still running, and on start/stop. That keeps
    Slurm ``.out`` files alive during long silent ``antechamber`` / ``g16``
    calls without dumping a table on every file refresh.
    """

    def __init__(
        self,
        store: JobProgressStore,
        *,
        board_path: Path,
        logger,
        interval_sec: float = 5.0,
        log_root_hint: str | None = None,
        thread_name: str = "job-progress-board",
        stream: TextIO | None = None,
        heartbeat_sec: float | None = None,
    ) -> None:
        import threading

        self.store = store
        self.board_path = Path(board_path)
        self.logger = logger
        self.interval_sec = float(interval_sec)
        self.log_root_hint = log_root_hint
        self.stream = stream
        if heartbeat_sec is None:
            self.heartbeat_sec = 20.0 if stream is not None else 0.0
        else:
            self.heartbeat_sec = float(heartbeat_sec)
        self._stop = threading.Event()
        self._thread = threading.Thread(
            target=self._loop, name=thread_name, daemon=True
        )
        self._last = ""
        self._last_fp = None
        self._last_stream_at = 0.0

    def start(self) -> None:
        self._emit(force_log=True)
        self._thread.start()

    def stop(self, *, final: bool = True) -> None:
        self._stop.set()
        self._thread.join(timeout=max(1.0, min(self.interval_sec, 2.0)))
        if final:
            self._emit(force_log=True)

    def _write_stream(self, text: str) -> None:
        out = self.stream
        if out is None:
            if self.logger is not None:
                self.logger.info("\n%s", text.rstrip("\n"))
            return
        try:
            payload = text if text.endswith("\n") else text + "\n"
            out.write("\n" + payload if not payload.startswith("\n") else payload)
            out.flush()
        except OSError:
            pass

    def _emit(self, *, force_log: bool = False) -> None:
        jobs = self.store.snapshot()
        text = format_job_board(
            jobs,
            title=self.store.title,
            id_header=self.store.id_header,
            empty_hint=self.store.empty_hint,
            detail_hint_label=self.store.detail_hint_label,
            log_root_hint=self.log_root_hint,
        )
        try:
            _atomic_write_text(self.board_path, text)
        except OSError:
            pass
        fp = jobs_fingerprint(jobs)
        now = time.monotonic()
        changed = fp != self._last_fp
        heartbeat_due = self.heartbeat_sec > 0 and (
            now - self._last_stream_at
        ) >= self.heartbeat_sec
        if force_log or changed or heartbeat_due:
            self._write_stream(text)
            self._last = text
            self._last_fp = fp
            self._last_stream_at = now

    def _loop(self) -> None:
        while not self._stop.wait(timeout=self.interval_sec):
            self._emit(force_log=False)


# ---------------------------------------------------------------------------
# Fragment / whole-ligand boards: one profile table, named aliases
# ---------------------------------------------------------------------------

from contextlib import contextmanager
from typing import Iterator

from .Console import attach_console_handlers, console_formatter, tee_stdio_to_file

_BOARD_PROFILES: dict[str, dict[str, Any]] = {
    "fragment": {
        "collection_key": "fragments",
        "id_header": "Fragment",
        "title": "Fragment dihedral twist - live status",
        "empty_hint": "no fragments registered yet",
        "detail_hint_label": "Per-fragment detail logs",
        "logger_ns": "ffpopt.workflows.frag",
        "watcher_thread": "frag-progress-board",
        "known_stages": (
            "queued",
            "prepare",
            "hl_scan",
            "orig_scan",
            "compare",
            "fit",
            "apply",
            "rescan",
            "finished",
            "failed",
        ),
    },
    "whole": {
        "collection_key": "batches",
        "id_header": "Batch",
        "title": "Whole-ligand dihedral twist - live status",
        "empty_hint": "no torsion batches registered yet",
        "detail_hint_label": "Per-batch detail logs",
        "logger_ns": "ffpopt.workflows.whole",
        "watcher_thread": "whole-progress-board",
        "known_stages": (
            "queued",
            "prepare",
            "twist",
            "finished",
            "failed",
        ),
    },
}


def _board_profile(kind: str) -> dict[str, Any]:
    try:
        return _BOARD_PROFILES[kind]
    except KeyError as exc:
        raise ValueError(f"unknown progress-board kind {kind!r}") from exc


def make_progress_store(kind: str, path: PathLike) -> JobProgressStore:
    """Build a :class:`JobProgressStore` from a named fragment/whole profile."""
    cfg = _board_profile(kind)
    return JobProgressStore(
        path,
        collection_key=cfg["collection_key"],
        id_header=cfg["id_header"],
        title=cfg["title"],
        empty_hint=cfg["empty_hint"],
        detail_hint_label=cfg["detail_hint_label"],
    )


def make_board_watcher(
    kind: str,
    store: JobProgressStore,
    *,
    board_path: Path,
    logger,
    interval_sec: float = 5.0,
    log_root_hint: str | None = None,
    stream: TextIO | None = None,
    heartbeat_sec: float | None = None,
) -> JobBoardWatcher:
    """Build a :class:`JobBoardWatcher` from a named fragment/whole profile."""
    cfg = _board_profile(kind)
    return JobBoardWatcher(
        store,
        board_path=board_path,
        logger=logger,
        interval_sec=interval_sec,
        log_root_hint=log_root_hint,
        thread_name=cfg["watcher_thread"],
        stream=stream,
        heartbeat_sec=heartbeat_sec,
    )


def format_kind_board(
    kind: str,
    jobs: Mapping[str, Mapping[str, Any]],
    *,
    title: str | None = None,
    log_root_hint: str | None = None,
) -> str:
    cfg = _board_profile(kind)
    return format_job_board(
        jobs,
        title=title if title is not None else cfg["title"],
        id_header=cfg["id_header"],
        empty_hint=cfg["empty_hint"],
        detail_hint_label=cfg["detail_hint_label"],
        log_root_hint=log_root_hint,
    )


def make_kind_file_logger(kind: str, job_id: str, log_path: Path):
    """Logger that writes to ``log_path`` and mirrors to the console."""
    import logging

    cfg = _board_profile(kind)
    name = f"{cfg['logger_ns']}.{job_id}"
    tag = f"ffpopt:{job_id}"
    logger = logging.getLogger(name)
    logger.handlers.clear()
    logger.propagate = False
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
    handler.setFormatter(console_formatter(tag))
    logger.addHandler(handler)
    attach_console_handlers(logger, tag=tag)
    return logger


@contextmanager
def kind_stdio_to_file(
    log_path: Path,
    *,
    job_id: str | None = None,
) -> Iterator[None]:
    """Tee stdout/stderr to ``log_path`` and the parent console."""
    tag = f"ffpopt:{job_id}" if job_id else "ffpopt"
    with tee_stdio_to_file(log_path, tag=tag):
        yield


KNOWN_FRAGMENT_STAGES = _BOARD_PROFILES["fragment"]["known_stages"]
KNOWN_STAGES = KNOWN_FRAGMENT_STAGES
KNOWN_WHOLE_STAGES = _BOARD_PROFILES["whole"]["known_stages"]


class FragmentProgressStore(JobProgressStore):
    """Process-shared fragment status table backed by a JSON file."""

    def __init__(self, path: PathLike) -> None:
        cfg = _board_profile("fragment")
        super().__init__(
            path,
            collection_key=cfg["collection_key"],
            id_header=cfg["id_header"],
            title=cfg["title"],
            empty_hint=cfg["empty_hint"],
            detail_hint_label=cfg["detail_hint_label"],
        )

    def register(
        self,
        fragment_id: str,
        *,
        bonds: int = 0,
        frag_dir: str | None = None,
        log_path: str | None = None,
    ) -> None:
        """Mark a fragment as queued before workers start."""
        detail = f"{bonds} bond(s)"
        if frag_dir:
            detail = f"{detail} | {frag_dir}"
        super().register(
            fragment_id,
            status="queued",
            stage="queued",
            detail=detail,
            bonds=bonds,
            log_path=log_path,
        )


def format_fragment_board(
    fragments: Mapping[str, Mapping[str, Any]],
    *,
    title: str = "Fragment dihedral twist - live status",
    log_root_hint: str | None = None,
) -> str:
    """Format ``{id: {status, stage, detail, ...}}`` as a fixed-width board."""
    return format_kind_board(
        "fragment", fragments, title=title, log_root_hint=log_root_hint
    )


@contextmanager
def fragment_stdio_to_file(
    log_path: Path,
    *,
    fragment_id: str | None = None,
) -> Iterator[None]:
    """Tee stdout/stderr to ``log_path`` and the parent console."""
    with kind_stdio_to_file(log_path, job_id=fragment_id):
        yield


def make_fragment_file_logger(fragment_id: str, log_path: Path):
    """Logger that writes to the fragment log and mirrors to the console."""
    return make_kind_file_logger("fragment", fragment_id, log_path)


class FragmentBoardWatcher(JobBoardWatcher):
    """Background thread that refreshes ``FRAG_STATUS.txt`` and logs on change."""

    def __init__(
        self,
        store: FragmentProgressStore,
        *,
        board_path: Path,
        logger,
        interval_sec: float = 5.0,
        log_root_hint: str | None = None,
    ) -> None:
        super().__init__(
            store,
            board_path=board_path,
            logger=logger,
            interval_sec=interval_sec,
            log_root_hint=log_root_hint,
            thread_name=_board_profile("fragment")["watcher_thread"],
        )


class WholeProgressStore(JobProgressStore):
    """Process-shared whole-ligand torsion-batch status table."""

    def __init__(self, path: PathLike) -> None:
        cfg = _board_profile("whole")
        super().__init__(
            path,
            collection_key=cfg["collection_key"],
            id_header=cfg["id_header"],
            title=cfg["title"],
            empty_hint=cfg["empty_hint"],
            detail_hint_label=cfg["detail_hint_label"],
        )

    def register(
        self,
        batch_id: str,
        *,
        bonds: int = 0,
        log_path: str | None = None,
    ) -> None:
        """Mark a torsion batch as queued before it starts."""
        super().register(
            batch_id,
            status="queued",
            stage="queued",
            detail=f"{bonds} bond(s)",
            bonds=bonds,
            log_path=log_path,
        )


def format_whole_board(
    batches: Mapping[str, Mapping[str, Any]],
    *,
    title: str = "Whole-ligand dihedral twist - live status",
    log_root_hint: str | None = None,
) -> str:
    """Format ``{id: {status, stage, detail, ...}}`` as a fixed-width board."""
    return format_kind_board(
        "whole", batches, title=title, log_root_hint=log_root_hint
    )


@contextmanager
def whole_stdio_to_file(
    log_path: Path,
    *,
    batch_id: str | None = None,
) -> Iterator[None]:
    """Tee stdout/stderr to ``log_path`` and the parent console."""
    with kind_stdio_to_file(log_path, job_id=batch_id):
        yield


def make_whole_file_logger(batch_id: str, log_path: Path):
    """Logger that writes to the batch log and mirrors to the console."""
    return make_kind_file_logger("whole", batch_id, log_path)


class WholeBoardWatcher(JobBoardWatcher):
    """Background thread that refreshes ``WHOLE_STATUS.txt`` and logs on change."""

    def __init__(
        self,
        store: WholeProgressStore,
        *,
        board_path: Path,
        logger,
        interval_sec: float = 5.0,
        log_root_hint: str | None = None,
    ) -> None:
        super().__init__(
            store,
            board_path=board_path,
            logger=logger,
            interval_sec=interval_sec,
            log_root_hint=log_root_hint,
            thread_name=_board_profile("whole")["watcher_thread"],
        )

