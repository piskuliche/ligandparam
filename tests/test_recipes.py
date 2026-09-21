"""Recipe registry and setup() stage graphs (no AmberTools / Gaussian)."""

from __future__ import annotations

import importlib
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

_TESTS_DIR = Path(__file__).resolve().parent
if str(_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(_TESTS_DIR))
import _paths  # noqa: E402


def setUpModule():
    _paths.ensure_ligandparam()


def _has_module(name: str) -> bool:
    try:
        importlib.import_module(name)
        return True
    except Exception:
        return False


def _require_rdkit(test: unittest.TestCase) -> None:
    if not _has_module("rdkit"):
        test.skipTest("rdkit required for recipe/stage imports")


class TestRecipeRegistry(unittest.TestCase):
    def test_available_recipes_matches_registry_keys(self):
        from ligandparam.recipes.Registry import _REGISTRY, available_recipes

        self.assertEqual(available_recipes(), sorted(_REGISTRY))
        for expected in (
            "freeligand",
            "lazyligand",
            "lazierligand",
            "dplazyligand",
            "dpfreeligand",
            "sqmligand",
        ):
            self.assertIn(expected, _REGISTRY)

    def test_unknown_recipe_raises(self):
        from ligandparam.recipes.Registry import get_recipe

        with self.assertRaises(ValueError) as ctx:
            get_recipe("not-a-recipe")
        self.assertIn("Unknown recipe", str(ctx.exception))


class TestRecipeSetupGraphs(unittest.TestCase):
    """Each registered recipe builds a non-empty, ordered stage list."""

    def setUp(self):
        _require_rdkit(self)

    def _tmp_recipe_args(self, td: str):
        cwd = Path(td)
        return cwd / "ligand.pdb", cwd

    def _assert_tail_parmchk_leap(self, stages):
        from ligandparam.stages import StageLeap, StageParmChk

        types = [type(s) for s in stages]
        self.assertEqual(types[-2], StageParmChk)
        self.assertEqual(types[-1], StageLeap)

    def test_freeligand_setup(self):
        from ligandparam.recipes.FreeLigand import FreeLigand
        from ligandparam.stages import (
            StageInitialize,
            StageMultiRespFit,
            StageNormalizeCharge,
        )

        with tempfile.TemporaryDirectory() as td:
            inp, cwd = self._tmp_recipe_args(td)
            recipe = FreeLigand(inp, cwd, net_charge=0, nproc=2, mem=4, logger="stream")
            recipe.setup()
            types = [type(s) for s in recipe.stages]
            self.assertGreaterEqual(len(types), 5)
            self.assertEqual(types[0], StageInitialize)
            self.assertFalse(recipe.stages[0].assign_charges)
            self.assertIn(StageNormalizeCharge, types)
            self.assertIn(StageMultiRespFit, types)
            self._assert_tail_parmchk_leap(recipe.stages)
            self.assertEqual(recipe.net_charge, 0)

    def test_freeligand_missing_net_charge(self):
        from ligandparam.recipes.FreeLigand import FreeLigand

        with self.assertRaises(KeyError):
            FreeLigand("ligand.pdb", "out_dir")

    def test_freeligand_bad_orientation_protocol(self):
        from ligandparam.recipes.FreeLigand import FreeLigand

        with self.assertRaises(ValueError):
            FreeLigand(
                "ligand.pdb",
                "out_dir",
                net_charge=0,
                orientation_protocol="not_a_protocol",
            )

    def test_lazyligand_setup(self):
        from ligandparam.recipes.LazyLigand import LazyLigand
        from ligandparam.stages import StageInitialize, StageLazyResp

        with tempfile.TemporaryDirectory() as td:
            inp, cwd = self._tmp_recipe_args(td)
            recipe = LazyLigand(inp, cwd, net_charge=-1, logger="stream")
            recipe.setup()
            types = [type(s) for s in recipe.stages]
            self.assertEqual(types[0], StageInitialize)
            self.assertIn(StageLazyResp, types)
            self._assert_tail_parmchk_leap(recipe.stages)

    def test_lazierligand_setup(self):
        from ligandparam.recipes.LazierLigand import LazierLigand
        from ligandparam.stages import StageInitialize, StageNormalizeCharge

        with tempfile.TemporaryDirectory() as td:
            inp, cwd = self._tmp_recipe_args(td)
            recipe = LazierLigand(inp, cwd, net_charge=0, logger="stream")
            recipe.setup()
            types = [type(s) for s in recipe.stages]
            self.assertEqual(types[0], StageInitialize)
            self.assertTrue(recipe.stages[0].assign_charges)
            self.assertIn(StageNormalizeCharge, types)
            self.assertEqual(sum(1 for t in types if t.__name__ == "StageParmChk"), 1)
            self._assert_tail_parmchk_leap(recipe.stages)

    def test_dpligand_setup_includes_dpminimize(self):
        from ligandparam.recipes.DpLazyLigand import DPLigand
        from ligandparam.stages import DPMinimize, StageInitialize

        with tempfile.TemporaryDirectory() as td:
            inp, cwd = self._tmp_recipe_args(td)
            recipe = DPLigand(inp, cwd, net_charge=0, logger="stream")
            recipe.setup()
            types = [type(s) for s in recipe.stages]
            self.assertEqual(types[0], StageInitialize)
            self.assertIn(DPMinimize, types)
            self._assert_tail_parmchk_leap(recipe.stages)

    def test_dpfreeligand_setup(self):
        from ligandparam.recipes.DpFreeLigand import DPFreeLigand
        from ligandparam.stages import DPMinimize, StageInitialize, StageMultiRespFit

        with tempfile.TemporaryDirectory() as td:
            inp, cwd = self._tmp_recipe_args(td)
            recipe = DPFreeLigand(inp, cwd, net_charge=0, logger="stream")
            recipe.setup()
            types = [type(s) for s in recipe.stages]
            self.assertEqual(types[0], StageInitialize)
            self.assertIn(DPMinimize, types)
            self.assertIn(StageMultiRespFit, types)
            self._assert_tail_parmchk_leap(recipe.stages)

    def test_sqmligand_setup(self):
        from ligandparam.recipes.OptLigand import SQMLigand
        from ligandparam.stages import StageInitialize, StageLazyResp
        from ligandparam.stages.DeepMd import DPMinimize

        with tempfile.TemporaryDirectory() as td:
            inp, cwd = self._tmp_recipe_args(td)
            recipe = SQMLigand(inp, cwd, net_charge=0, logger="stream")
            recipe.setup()
            types = [type(s) for s in recipe.stages]
            self.assertEqual(types[0], StageInitialize)
            self.assertIn(StageLazyResp, types)
            self.assertNotIn(DPMinimize, types)
            self._assert_tail_parmchk_leap(recipe.stages)

    def test_lazierligand_execute_forwards_overrides(self):
        from ligandparam.recipes.LazierLigand import LazierLigand

        with tempfile.TemporaryDirectory() as td:
            inp, cwd = self._tmp_recipe_args(td)
            recipe = LazierLigand(inp, cwd, net_charge=0, nproc=4, logger="stream")
            recipe.setup()
            seen = []

            def _capture(stage):
                def _exec(*, dry_run=False, nproc=None, mem=None):
                    seen.append((dry_run, nproc, mem))

                return _exec

            for stage in recipe.stages:
                stage.execute = _capture(stage)
            recipe.execute(dry_run=True, nproc=8, mem=16)
            self.assertTrue(seen)
            self.assertTrue(all(t == (True, 8, 16) for t in seen))

    def test_dihed_correct_does_not_append_twist_stage(self):
        """Recipes record dihed_correct; ALPS runs the twist."""
        from ligandparam.recipes.FreeLigand import FreeLigand

        with tempfile.TemporaryDirectory() as td:
            inp, cwd = self._tmp_recipe_args(td)
            recipe = FreeLigand(
                inp,
                cwd,
                net_charge=0,
                dihed_correct=True,
                dihed_model="xtb",
                dihed_delta=15,
                logger="stream",
            )
            recipe.setup()
            self.assertTrue(recipe.dihed_correct)
            self.assertEqual(recipe.dihed_model, "xtb")
            self.assertEqual(recipe.dihed_delta, 15)
            names = [type(s).__name__ for s in recipe.stages]
            self.assertNotIn("StageDihedTwistCorrection", names)

    def test_dry_run_execute_invokes_each_stage(self):
        from ligandparam.recipes.LazierLigand import LazierLigand

        with tempfile.TemporaryDirectory() as td:
            inp, cwd = self._tmp_recipe_args(td)
            recipe = LazierLigand(inp, cwd, net_charge=0, logger="stream")
            recipe.setup()
            calls = []
            for stage in recipe.stages:
                stage.execute = MagicMock(
                    side_effect=lambda *a, _s=stage, **k: calls.append(_s.stage_name)
                )
            recipe.execute(dry_run=True)
            self.assertEqual(len(calls), len(recipe.stages))
            self.assertEqual(calls, [s.stage_name for s in recipe.stages])

    def test_every_registry_recipe_uses_common_builders(self):
        from ligandparam.recipes.Registry import _REGISTRY, get_recipe

        for name, path_cls in _REGISTRY.items():
            mod_path = path_cls.split(":")[0]
            mod = importlib.import_module(mod_path)
            src = Path(mod.__file__).read_text(encoding="utf-8")
            self.assertIn(
                "ligandparam.recipes.Common",
                src,
                f"{name} should import recipes.common builders",
            )
            with tempfile.TemporaryDirectory() as td:
                inp, cwd = self._tmp_recipe_args(td)
                recipe = get_recipe(
                    name,
                    in_filename=str(inp),
                    cwd=str(cwd),
                    net_charge=0,
                    logger="stream",
                )
                recipe.setup()
                self.assertGreater(len(recipe.stages), 0, f"{name} setup() empty")
                self._assert_tail_parmchk_leap(recipe.stages)


class TestCommonRecipeTail(unittest.TestCase):
    def test_charge_update_parmchk_leap_order(self):
        _require_rdkit(self)
        from ligandparam.recipes.Common import charge_update_parmchk_leap_stages
        from ligandparam.stages import StageLeap, StageParmChk, StageUpdate

        recipe = SimpleNamespace(cwd=Path("."), net_charge=0, logger=None, kwargs={})
        stages = charge_update_parmchk_leap_stages(
            recipe=recipe,
            initial_mol2="a.mol2",
            final_mol2="b.mol2",
            nonminimized_mol2="c.mol2",
            frcmod="x.frcmod",
            lib="x.lib",
        )
        self.assertEqual(
            [type(s) for s in stages],
            [StageUpdate, StageParmChk, StageLeap],
        )

    def test_init_normalize_center_and_gaussian_kwargs(self):
        _require_rdkit(self)
        from ligandparam.recipes.Common import (
            gaussian_runtime_kwargs,
            init_normalize_center_stages,
            rotation_stage_kwargs,
        )
        from ligandparam.stages import StageDisplaceMol, StageInitialize, StageNormalizeCharge

        recipe = SimpleNamespace(
            in_filename=Path("lig.pdb"),
            cwd=Path("."),
            net_charge=0,
            logger=None,
            kwargs={},
            nproc=2,
            mem=4,
            gaussian_root=None,
            gauss_exedir=None,
            gaussian_binary=None,
            gaussian_scratch=None,
            force_gaussian_rerun=False,
            orientation_protocol="so3_n28",
            theory={"low": "HF/6-31G*", "high": "PBE1PBE/6-31G*"},
        )
        stages = init_normalize_center_stages(
            recipe=recipe,
            initial_mol2="i.mol2",
            centered_out="c.mol2",
        )
        self.assertEqual(
            [type(s) for s in stages],
            [StageInitialize, StageNormalizeCharge, StageDisplaceMol],
        )
        self.assertFalse(stages[0].assign_charges)
        gkw = gaussian_runtime_kwargs(recipe)
        self.assertEqual(gkw["nproc"], 2)
        self.assertEqual(gkw["mem"], 4)
        self.assertIn("force_gaussian_rerun", gkw)
        rkw = rotation_stage_kwargs(recipe)
        self.assertEqual(rkw["orientation_protocol"], "so3_n28")


class TestDihedOptions(unittest.TestCase):
    def test_pop_and_apply_dihed_options(self):
        from ligandparam.recipes.DihedOptions import apply_dihed_options, pop_dihed_options

        kwargs = {
            "dihed_correct": True,
            "dihed_model": "xtb",
            "dihed_delta": 5,
            "keep": 1,
        }
        opts = pop_dihed_options(dict(kwargs))
        self.assertTrue(opts["dihed_correct"])
        self.assertEqual(opts["dihed_delta"], 5)
        obj = SimpleNamespace()
        apply_dihed_options(obj, kwargs)
        self.assertTrue(obj.dihed_correct)
        self.assertEqual(obj.dihed_model, "xtb")
        self.assertEqual(kwargs, {"keep": 1})

    def test_coerce_fragment_config(self):
        from ligandparam.recipes.DihedOptions import coerce_fragment_config

        self.assertIsNone(coerce_fragment_config(None))
        cfg = coerce_fragment_config({"angle_step": 15, "cap_strategy": "hydrogen"})
        self.assertEqual(cfg["angle_step"], 15)
        self.assertEqual(coerce_fragment_config("pass"), "pass")

    def test_dihed_fragment_strategy_merges_into_config(self):
        from ligandparam.recipes.DihedOptions import pop_dihed_options

        opts = pop_dihed_options({"dihed_fragment_strategy": "pfizer"})
        self.assertEqual(opts["dihed_fragment_strategy"], "pfizer")
        self.assertEqual(opts["dihed_fragment_config"]["strategy"], "pfizer")
        opts2 = pop_dihed_options(
            {
                "dihed_fragment_config": {"angle_step": 15},
                "dihed_fragment_strategy": "wbo",
            }
        )
        self.assertEqual(opts2["dihed_fragment_config"]["strategy"], "wbo")
        self.assertEqual(opts2["dihed_fragment_config"]["angle_step"], 15)


class TestRecipeDefaultsIsolation(unittest.TestCase):
    def test_fresh_defaults_not_shared(self):
        from ligandparam.Parametrization import fresh_recipe_defaults

        a = fresh_recipe_defaults()
        b = fresh_recipe_defaults()
        a["leaprc"].append("leaprc.protein.ff14SB")
        self.assertNotIn("leaprc.protein.ff14SB", b["leaprc"])
        a["theory"]["low"] = "X"
        self.assertNotEqual(b["theory"]["low"], "X")


if __name__ == "__main__":
    unittest.main()
