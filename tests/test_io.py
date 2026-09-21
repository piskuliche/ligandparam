"""Amber bundle resolution and logging contracts."""

from __future__ import annotations

import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

_TESTS_DIR = Path(__file__).resolve().parent
if str(_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(_TESTS_DIR))
import _paths  # noqa: E402


def setUpModule():
    _paths.ensure_ligandparam()


class TestAmberBundleIO(unittest.TestCase):
    def _touch_triplet(self, work_dir: Path, stem: str) -> None:
        (work_dir / f"{stem}.mol2").write_text("@<TRIPOS>MOLECULE\n", encoding="utf-8")
        (work_dir / f"{stem}.lib").write_text("!entry\n", encoding="utf-8")
        (work_dir / f"{stem}.frcmod").write_text("Remark line\n", encoding="utf-8")

    def test_resolve_explicit_paths(self):
        from ligandparam.io.AmberBundle import AmberLigandBundle, resolve_getparam_bundle

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            self._touch_triplet(root, "chaps")
            bundle = resolve_getparam_bundle(
                mol2=root / "chaps.mol2",
                lib=root / "chaps.lib",
                frcmod=root / "chaps.frcmod",
            )
            self.assertIsInstance(bundle, AmberLigandBundle)
            self.assertEqual(bundle.stem, "chaps")
            self.assertEqual(bundle.work_dir, root.resolve())

    def test_resolve_getparam_layout_with_label(self):
        from ligandparam.io.AmberBundle import resolve_getparam_bundle

        with tempfile.TemporaryDirectory() as td:
            cwd = Path(td)
            work = cwd / "CHA3" / "CHA"
            work.mkdir(parents=True)
            self._touch_triplet(work, "chaps")
            bundle = resolve_getparam_bundle(
                cwd=cwd, data_cwd="CHA3", resname="CHA", label="chaps"
            )
            self.assertEqual(bundle.stem, "chaps")
            self.assertEqual(bundle.work_dir, work.resolve())

    def test_resolve_label_is_case_insensitive(self):
        from ligandparam.io.AmberBundle import resolve_getparam_bundle

        with tempfile.TemporaryDirectory() as td:
            cwd = Path(td)
            work = cwd / "SDS" / "SDS"
            work.mkdir(parents=True)
            self._touch_triplet(work, "SDS")
            bundle = resolve_getparam_bundle(
                cwd=cwd, data_cwd="SDS", resname="SDS", label="sds"
            )
            self.assertEqual(bundle.mol2.name, "SDS.mol2")
            self.assertEqual(bundle.lib.name, "SDS.lib")
            self.assertEqual(bundle.frcmod.name, "SDS.frcmod")

    def test_resolve_unique_triplet_if_label_missing(self):
        from ligandparam.io.AmberBundle import resolve_getparam_bundle

        with tempfile.TemporaryDirectory() as td:
            cwd = Path(td)
            work = cwd / "CHA3" / "CHA"
            work.mkdir(parents=True)
            self._touch_triplet(work, "chaps")
            (work / "chaps.initial.mol2").write_text("@<TRIPOS>MOLECULE\n", encoding="utf-8")
            bundle = resolve_getparam_bundle(
                cwd=cwd, data_cwd="CHA3", resname="CHA"
            )
            self.assertEqual(bundle.stem, "chaps")

    def test_missing_triplet_raises(self):
        from ligandparam.io.AmberBundle import resolve_getparam_bundle

        with tempfile.TemporaryDirectory() as td:
            cwd = Path(td)
            (cwd / "A" / "B").mkdir(parents=True)
            with self.assertRaises(FileNotFoundError):
                resolve_getparam_bundle(cwd=cwd, data_cwd="A", resname="B", label="x")

    def test_to_scission_input(self):
        from ligandparam.io.AmberBundle import resolve_getparam_bundle

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            self._touch_triplet(root, "LIG")
            bundle = resolve_getparam_bundle(
                mol2=root / "LIG.mol2",
                lib=root / "LIG.lib",
                frcmod=root / "LIG.frcmod",
            )
            inp = bundle.to_scission_input()
            self.assertEqual(inp["mol2_path"], bundle.mol2)
            self.assertEqual(inp["lib_path"], bundle.lib)
            self.assertEqual(inp["frcmod_path"], bundle.frcmod)


class TestLogging(unittest.TestCase):
    def test_success_quote_skips_comments_and_empty(self):
        from ligandparam.Log import (
            dihed_correct_ok,
            format_reminder_line,
            load_quotes,
            log_success_quote,
        )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "quotes.txt"
            path.write_text(
                "# comment\n\n  first quote  \n\"wrapped\"\n",
                encoding="utf-8",
            )
            self.assertEqual(load_quotes(path), ["first quote", "wrapped"])
            empty = root / "empty.txt"
            empty.write_text("# only comments\n\n", encoding="utf-8")
            self.assertEqual(load_quotes(empty), [])
            self.assertIsNone(log_success_quote(quotes_path=empty))

            with patch.dict("os.environ", {}, clear=False):
                import os

                os.environ.pop("ALPS_BANNER_PRINTED", None)
                self.assertEqual(
                    format_reminder_line("hello"),
                    "LIGANDPARAM reminds you: hello",
                )
                self.assertEqual(
                    format_reminder_line(
                        '"It is nothing to die; it is dreadful not to live."'
                        " - Victor Hugo, Les Miserables"
                    ),
                    'LIGANDPARAM reminds you: "It is nothing to die; '
                    'it is dreadful not to live." - Victor Hugo, Les Miserables',
                )
                with patch.dict("os.environ", {"ALPS_BANNER_PRINTED": "1"}):
                    self.assertEqual(
                        format_reminder_line("hello"),
                        "ALPS reminds you: hello",
                    )
            self.assertEqual(
                format_reminder_line("hello", speaker="ALPS"),
                "ALPS reminds you: hello",
            )
            self.assertEqual(
                format_reminder_line("hello", speaker="LIGANDPARAM"),
                "LIGANDPARAM reminds you: hello",
            )
            self.assertFalse(dihed_correct_ok(None))
            self.assertFalse(dihed_correct_ok({"merged_frcmod": "/no/such.frcmod"}))
            frc = root / "out.frcmod"
            frc.write_text("DIHE\n", encoding="utf-8")
            self.assertTrue(dihed_correct_ok({"merged_frcmod": str(frc)}))
            self.assertFalse(
                dihed_correct_ok(
                    {
                        "merged_frcmod": str(frc),
                        "fragments": [{"fragment_id": "f1", "status": "failed"}],
                    }
                )
            )
            self.assertFalse(dihed_correct_ok({"merged_frcmod": str(frc)}, dry_run=True))
            with patch.dict("os.environ", {}, clear=False):
                import os

                os.environ.pop("ALPS_BANNER_PRINTED", None)
                with patch("sys.stdout", io.StringIO()) as buf:
                    picked = log_success_quote(quotes_path=path)
            self.assertIn(picked, ("first quote", "wrapped"))
            out = buf.getvalue()
            self.assertIn("LIGANDPARAM reminds you:", out)
            self.assertNotIn('reminds you: "', out)


if __name__ == "__main__":
    unittest.main()
