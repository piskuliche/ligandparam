"""Antechamber argv must not treat a negative net charge as another flag."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock

_TESTS_DIR = Path(__file__).resolve().parent
if str(_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(_TESTS_DIR))
import _paths  # noqa: E402


def setUpModule():
    _paths.ensure_ligandparam()


class TestAntechamberNetChargeArgv(unittest.TestCase):
    def test_negative_nc_is_not_a_bare_minus_token(self):
        from ligandparam.Interfaces import _format_program_flag

        flag, value = _format_program_flag("nc", -1.0)
        self.assertEqual(flag, "-nc")
        self.assertFalse(value.startswith("-"), value)
        self.assertEqual(int(value), -1)

    def test_zero_and_positive_nc_stay_plain_ints(self):
        from ligandparam.Interfaces import _format_program_flag

        self.assertEqual(_format_program_flag("nc", 0.0), ["-nc", "0"])
        self.assertEqual(_format_program_flag("nc", 2), ["-nc", "2"])

    def test_antechamber_dry_run_logs_spaced_negative_nc(self):
        from ligandparam.Interfaces import Antechamber

        logger = MagicMock()
        ante = Antechamber(cwd=".", logger=logger, nproc=1)
        ante.call(i="lig.mol2", fi="mol2", o="out.mol2", fo="mol2", c="bcc", nc=-1.0, dry_run=True)
        logged = " ".join(str(c) for call in logger.info.call_args_list for c in call.args)
        self.assertIn("-nc", logged)
        self.assertNotIn("-nc -1", logged)
        self.assertIn("-1", logged)


class TestGaussianFailureIncludesLog(unittest.TestCase):
    def test_failure_message_appends_log_tail(self):
        from types import SimpleNamespace

        from ligandparam.Interfaces import _subprocess_failure_message

        with tempfile.TemporaryDirectory() as td:
            cwd = Path(td)
            log = cwd / "job.log"
            log.write_text(
                "Entering Link 1\nError termination via Lnk1e\n"
                "The combination of multiplicity 1 and 109 electrons is impossible.\n",
                encoding="utf-8",
            )
            p = SimpleNamespace(returncode=1, stdout=b"", stderr=b"")
            msg = _subprocess_failure_message(cwd, p, tool="Gaussian", log_path=log)
            self.assertIn("returncode=1", msg)
            self.assertIn("109 electrons is impossible", msg)
