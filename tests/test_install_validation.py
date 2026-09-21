"""Install validation for standalone ligandparam.

Run after ``pip install -e .`` (or from this checkout)::

    python -m unittest tests.test_install_validation -v
"""

from __future__ import annotations

import importlib
import importlib.metadata
import importlib.util
import sys
import unittest
from pathlib import Path

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


class TestCorePackageInstall(unittest.TestCase):
    def test_version(self):
        import ligandparam

        self.assertTrue(ligandparam.__version__)
        self.assertRegex(ligandparam.__version__, r"^\d+\.\d+")
        self.assertEqual(ligandparam.__logging_name__, "ligandparam")


class TestPublicAPISurface(unittest.TestCase):
    def test_stages_and_recipes(self):
        from ligandparam.recipes.Registry import available_recipes, get_recipe
        from ligandparam.stages.AbstractStage import AbstractStage

        names = available_recipes()
        self.assertIsInstance(names, (list, tuple))
        self.assertIn("freeligand", names)
        self.assertIn("lazyligand", names)
        with self.assertRaises(ValueError):
            get_recipe("not-a-real-recipe")
        self.assertTrue(issubclass(AbstractStage, object))

        for modname in (
            "ligandparam.recipes.FreeLigand",
            "ligandparam.recipes.LazyLigand",
            "ligandparam.stages.Gaussian",
        ):
            self.assertIsNotNone(importlib.util.find_spec(modname), modname)

        m = importlib.import_module("ligandparam.stages.Gaussian")
        for name in (
            "GaussianMinimizeRESP",
            "GaussianRESP",
            "StageGaussianRotation",
            "StageGaussianToMol2",
            "StageGaussiantoMol2",
        ):
            self.assertTrue(isinstance(getattr(m, name), type), name)
        self.assertIs(m.StageGaussiantoMol2, m.StageGaussianToMol2)

        if _has_module("rdkit"):
            pdb_names = importlib.import_module("ligandparam.stages.PdbNames")
            self.assertIs(pdb_names.PDB_Name_Fixer, pdb_names.StagePdbNameFixer)

    def test_smiles_stage(self):
        if not _has_module("rdkit"):
            self.skipTest("rdkit required for smiles stages")
        smiles = importlib.import_module("ligandparam.stages.SmilesToPdb")
        self.assertTrue(hasattr(smiles, "StageSmilesToPDB"))


class TestCLIEntrypoints(unittest.TestCase):
    EXPECTED = (
        ("ligandparam.cli.LigGetParam", "main"),
        ("ligandparam.cli.SmilesToPdb", "main"),
        ("ligandparam.cli.CliLigHFix", "lighfix"),
    )

    def test_cli_modules_expose_callables(self):
        if not _has_module("rdkit"):
            self.skipTest("rdkit required for ligandparam CLIs")
        for modname, attr in self.EXPECTED:
            with self.subTest(module=modname):
                mod = importlib.import_module(modname)
                self.assertTrue(callable(getattr(mod, attr)))

    def test_installed_console_scripts_when_available(self):
        try:
            dist = importlib.metadata.distribution("ligandparam")
        except importlib.metadata.PackageNotFoundError:
            self.skipTest("ligandparam distribution metadata unavailable")
        ep_names = {ep.name for ep in dist.entry_points if ep.group == "console_scripts"}
        for required in ("lig-getparam", "lighfix", "smiles-to-pdb", "lig-to-sage"):
            self.assertIn(required, ep_names, required)
        for banned in ("lig-dihed-correct", "lig-scission"):
            self.assertNotIn(banned, ep_names, banned)


class TestIsolation(unittest.TestCase):
    def test_source_does_not_import_alps_ffpopt_or_scission(self):
        import ast

        root = _paths.package_root()
        hits: list[str] = []
        banned = frozenset({"alps", "ffpopt", "scission"})
        for fp in root.rglob("*.py"):
            if "__pycache__" in fp.parts or "tests" in fp.parts:
                continue
            text = fp.read_text(encoding="utf-8", errors="replace")
            try:
                tree = ast.parse(text, filename=str(fp))
            except SyntaxError:
                continue
            rel = fp.relative_to(root)
            for node in ast.walk(tree):
                names: list[str] = []
                if isinstance(node, ast.Import):
                    names = [alias.name for alias in node.names]
                elif isinstance(node, ast.ImportFrom) and node.module:
                    names = [node.module]
                for name in names:
                    if name.split(".")[0] in banned:
                        hits.append(f"{rel}:{node.lineno} {name}")
        self.assertEqual(hits, [], "ligandparam imports companions:\n" + "\n".join(hits))


if __name__ == "__main__":
    unittest.main()
