"""Stage contracts that do not need AmberTools or Gaussian."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

_TESTS_DIR = Path(__file__).resolve().parent
if str(_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(_TESTS_DIR))
import _paths  # noqa: E402


def setUpModule():
    _paths.ensure_ligandparam()


class TestStageChargeNormalize(unittest.TestCase):
    def _stage(self, net_charge=0, precision=0.001, decimals=3):
        from ligandparam.stages.Charge import StageNormalizeCharge

        st = StageNormalizeCharge.__new__(StageNormalizeCharge)
        st.net_charge = net_charge
        st.precision = precision
        st.decimals = decimals
        st.logger = MagicMock()
        return st

    def test_nonzero_net_charge(self):
        st = self._stage(net_charge=1, precision=0.001, decimals=3)
        q = [0.4, 0.3, 0.2]
        rounded, total, diff = st.check_charge(q)
        out = st.normalize(rounded, diff)
        _, new_total, _ = st.check_charge(out)
        self.assertTrue(abs(new_total - 1.0) < 0.002)

    def test_zero_count_safe(self):
        st = self._stage(net_charge=0)
        out = st.normalize([0.0, 0.0], 0.0)
        self.assertEqual(list(out), [0.0, 0.0])

    def test_large_delta_warns(self):
        st = self._stage(net_charge=1)
        out = st.normalize([0.0, 0.0], 0.05)
        self.assertAlmostEqual(float(sum(out)), 0.05, places=6)
        st.logger.warning.assert_called()


class TestAbstractStageTemplate(unittest.TestCase):
    def test_execute_calls_run_and_tracks_new_files(self):
        from ligandparam.stages.AbstractStage import AbstractStage

        class _Tiny(AbstractStage):
            def _run(self, dry_run=False, nproc=None, mem=None):
                self.seen = (dry_run, nproc, mem)
                (self.cwd / "created.txt").write_text("x", encoding="utf-8")
                return "ok"

        with tempfile.TemporaryDirectory() as td:
            cwd = Path(td)
            stage = _Tiny("Tiny", cwd / "in.pdb", cwd, logger=MagicMock())
            stage.required = []
            out = stage.execute(dry_run=True, nproc=3, mem=7)
            self.assertEqual(out, "ok")
            self.assertEqual(stage.seen, (True, 3, 7))
            self.assertIn("created.txt", stage.new_files)


class TestCpuBudget(unittest.TestCase):
    def test_gaussian_orientation_budget_splits_memory(self):
        from ligandparam.runtime.CpuBudget import split_gaussian_orientation_budget

        n_workers, job_nproc, job_mem = split_gaussian_orientation_budget(60, 28, 32)
        self.assertGreaterEqual(job_mem, 4)
        self.assertLessEqual(n_workers * job_nproc, 60)
        self.assertLessEqual(n_workers * job_mem, 32)


class TestGaussianInterface(unittest.TestCase):
    def test_gaussian_failure_includes_returncode(self):
        from ligandparam.Interfaces import Gaussian

        class _Proc:
            returncode = 137
            stdout = b""
            stderr = b"Killed"

        with tempfile.TemporaryDirectory() as tmp:
            gau = Gaussian(
                cwd=tmp,
                gaussian_root="",
                gauss_exedir="",
                gaussian_binary="g16",
                gaussian_scratch="",
            )
            with patch("ligandparam.Interfaces.subprocess.run", return_value=_Proc()):
                with self.assertRaises(RuntimeError) as ctx:
                    gau.call(inp_pipe="job.com", out_pipe="job.log")
        msg = str(ctx.exception)
        self.assertIn("returncode=137", msg)
        self.assertIn("OOM", msg)


class TestParmHelperExports(unittest.TestCase):
    def test_parmhelper_saveparm_exported(self):
        import ast

        root = _paths.package_root()
        utils = root / "multiresp" / "ParmEdUtils.py"
        helper = root / "multiresp" / "ParmHelper.py"
        util_names = {
            n.name
            for n in ast.parse(utils.read_text(encoding="utf-8")).body
            if isinstance(n, ast.FunctionDef)
        }
        self.assertIn("SaveParm", util_names)
        text = helper.read_text(encoding="utf-8")
        self.assertIn("SaveParm", text)
        self.assertIn("from ligandparam.multiresp.ParmEdUtils import", text)

    def test_parmhelper_sed_escape_is_valid(self):
        import warnings

        src = _paths.package_root() / "multiresp" / "ParmHelper.py"
        text = src.read_text(encoding="utf-8")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", SyntaxWarning)
            compile(text, str(src), "exec")
        syn = [w for w in caught if issubclass(w.category, SyntaxWarning)]
        self.assertEqual(syn, [], [str(w.message) for w in syn])


if __name__ == "__main__":
    unittest.main()
