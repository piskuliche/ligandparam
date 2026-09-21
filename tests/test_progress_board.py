"""Live recipe / job-status board (no AmberTools / Gaussian)."""

from __future__ import annotations

import io
import logging
import sys
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace

_TESTS_DIR = Path(__file__).resolve().parent
if str(_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(_TESTS_DIR))
import _paths  # noqa: E402


def setUpModule():
    _paths.ensure_ligandparam()


class TestBoardFormatting(unittest.TestCase):
    def test_elapsed_and_fingerprint_ignore_clock(self):
        from ligandparam.runtime.ProgressBoard import (
            elapsed_phrase,
            format_job_board,
            jobs_fingerprint,
        )

        started = time.strftime(
            "%Y-%m-%d %H:%M:%S", time.localtime(time.time() - 95)
        )
        jobs = {
            "Initialize": {
                "status": "running",
                "stage": "running",
                "detail": "StageInitialize",
                "started": started,
            }
        }
        board = format_job_board(jobs, title="recipe live")
        self.assertIn("recipe live", board)
        self.assertIn("refreshed", board)
        self.assertIn("elapsed", board)
        self.assertRegex(elapsed_phrase(started), r"1m")
        fp = jobs_fingerprint(jobs)
        jobs["Initialize"]["started"] = time.strftime("%Y-%m-%d %H:%M:%S")
        self.assertEqual(fp, jobs_fingerprint(jobs))

    def test_running_sets_started_once(self):
        from ligandparam.runtime.ProgressBoard import JobProgressStore

        with tempfile.TemporaryDirectory() as td:
            store = JobProgressStore(Path(td) / "p.json")
            store.register("A", status="queued", stage="queued")
            store.update("A", status="running", stage="running")
            first = store.snapshot()["A"]["started"]
            first_epoch = store.snapshot()["A"]["started_epoch"]
            store.update("A", status="running", stage="running")
            snap = store.snapshot()["A"]
            self.assertEqual(snap["started"], first)
            self.assertEqual(snap["started_epoch"], first_epoch)

    def test_new_stage_restarts_elapsed_clock(self):
        from ligandparam.runtime.ProgressBoard import JobProgressStore

        with tempfile.TemporaryDirectory() as td:
            store = JobProgressStore(Path(td) / "p.json")
            store.update("A", status="running", stage="prepare")
            first = store.snapshot()["A"]["started_epoch"]
            time.sleep(0.05)
            store.update("A", status="running", stage="twist")
            self.assertGreater(store.snapshot()["A"]["started_epoch"], first)

    def test_reregister_queued_clears_stale_started(self):
        from ligandparam.runtime.ProgressBoard import JobProgressStore

        with tempfile.TemporaryDirectory() as td:
            store = JobProgressStore(Path(td) / "p.json")
            store.update("A", status="running", stage="twist")
            self.assertIn("started_epoch", store.snapshot()["A"])
            store.register("A", status="queued", stage="queued")
            snap = store.snapshot()["A"]
            self.assertNotIn("started", snap)
            self.assertNotIn("started_epoch", snap)

    def test_elapsed_prefers_epoch(self):
        from ligandparam.runtime.ProgressBoard import elapsed_phrase, format_job_board

        stale = time.strftime(
            "%Y-%m-%d %H:%M:%S", time.localtime(time.time() - 273)
        )
        self.assertRegex(elapsed_phrase(stale, time.time()), r"^\d+s$")
        board = format_job_board(
            {
                "fragment_1": {
                    "status": "running",
                    "stage": "twist",
                    "detail": "2 bond(s)",
                    "started": stale,
                    "started_epoch": time.time(),
                }
            }
        )
        self.assertIn("elapsed", board)
        self.assertNotIn("4m", board)

    def test_fragments_list_in_numeric_order(self):
        from ligandparam.runtime.ProgressBoard import format_fragment_board

        board = format_fragment_board(
            {
                "fragment_10": {"status": "running", "stage": "twist", "detail": "1 bond(s)"},
                "fragment_2": {"status": "done", "stage": "finished", "detail": "1 bond(s)"},
                "fragment_1": {"status": "done", "stage": "finished", "detail": "2 bond(s)"},
                "fragment_11": {"status": "queued", "stage": "queued", "detail": "1 bond(s)"},
            }
        )
        names = [
            line.split()[0]
            for line in board.splitlines()
            if line.strip().startswith("fragment_")
        ]
        self.assertEqual(
            names,
            ["fragment_1", "fragment_2", "fragment_10", "fragment_11"],
        )


class TestJobBoardWatcher(unittest.TestCase):
    def test_stream_heartbeat_and_status_file(self):
        from ligandparam.runtime.ProgressBoard import JobBoardWatcher, JobProgressStore

        with tempfile.TemporaryDirectory() as td:
            cwd = Path(td)
            store = JobProgressStore(cwd / "p.json", title="Heartbeat board")
            store.register(
                "Initialize",
                status="running",
                stage="running",
                detail="StageInitialize",
            )
            buf = io.StringIO()
            logger = logging.getLogger("ligandparam.test.board")
            logger.handlers.clear()
            logger.addHandler(logging.NullHandler())
            watcher = JobBoardWatcher(
                store,
                board_path=cwd / "STATUS.txt",
                logger=logger,
                interval_sec=30.0,
                heartbeat_sec=0.01,
                stream=buf,
            )
            watcher.start()
            try:
                n0 = buf.getvalue().count("Heartbeat board")
                self.assertGreaterEqual(n0, 1)
                watcher._last_stream_at = 0.0
                watcher._emit(force_log=False)
                self.assertGreater(buf.getvalue().count("Heartbeat board"), n0)
            finally:
                watcher.stop()
            text = buf.getvalue()
            self.assertIn("Initialize", text)
            self.assertTrue((cwd / "STATUS.txt").is_file())
            self.assertIn("Heartbeat board", (cwd / "STATUS.txt").read_text(encoding="utf-8"))


class TestDriverRecipeBoard(unittest.TestCase):
    def test_execute_writes_recipe_status_and_marks_done(self):
        from ligandparam.Driver import Driver

        class _Drv(Driver):
            def __init__(self, cwd):
                self.cwd = Path(cwd)
                self.logger = logging.getLogger("ligandparam.test.driver")
                self.progress_stream = io.StringIO()
                self.stages = [
                    SimpleNamespace(
                        stage_name="Tiny",
                        execute=lambda **k: None,
                    ),
                    SimpleNamespace(
                        stage_name="Tiny",
                        execute=lambda **k: None,
                    ),
                ]

        with tempfile.TemporaryDirectory() as td:
            drv = _Drv(td)
            drv.execute(dry_run=True)
            board = Path(td) / "RECIPE_STATUS.txt"
            self.assertTrue(board.is_file())
            text = board.read_text(encoding="utf-8")
            self.assertIn("Tiny", text)
            self.assertIn("Tiny#2", text)
            self.assertIn("done", text)
            streamed = drv.progress_stream.getvalue()
            self.assertIn("ligandparam recipe", streamed)
            snap = (Path(td) / ".recipe_progress.json").read_text(encoding="utf-8")
            self.assertIn('"status": "done"', snap)

    def test_failed_stage_is_marked_failed(self):
        from ligandparam.Driver import Driver

        class _Drv(Driver):
            def __init__(self, cwd):
                self.cwd = Path(cwd)
                self.logger = logging.getLogger("ligandparam.test.driver")
                self.progress_stream = io.StringIO()
                self.stages = [
                    SimpleNamespace(
                        stage_name="Boom",
                        execute=_raise,
                    )
                ]

        with tempfile.TemporaryDirectory() as td:
            drv = _Drv(td)
            with self.assertRaises(RuntimeError):
                drv.execute()
            snap = (Path(td) / ".recipe_progress.json").read_text(encoding="utf-8")
            self.assertIn('"status": "failed"', snap)


def _raise(**_k):
    raise ValueError("nope")


if __name__ == "__main__":
    unittest.main()
