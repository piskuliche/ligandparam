"""Import-level smoke tests.

These need no external tools (no Gaussian, no AmberTools) and are meant to catch the
class of defect that shipped in 1.0.0: a console script pointing at a function that
does not exist, and a runtime dependency that was only ever satisfied transitively.
"""

import importlib
import importlib.metadata
import pkgutil

import pytest

import ligandparam


def _module_names():
    """Every importable module in the package, excluding the unmaintained ones."""
    names = []
    for info in pkgutil.walk_packages(ligandparam.__path__, "ligandparam."):
        # `deprecated/` imports pre-reorganisation module paths and is excluded from
        # the wheel; it is kept in the repo for reference only.
        if ".deprecated" in info.name:
            continue
        names.append(info.name)
    return sorted(names)


MODULE_NAMES = _module_names()


def test_found_modules():
    assert len(MODULE_NAMES) > 20, f"walk_packages only found {MODULE_NAMES}"


@pytest.mark.parametrize("module_name", MODULE_NAMES)
def test_module_imports(module_name):
    importlib.import_module(module_name)


def _console_scripts():
    eps = importlib.metadata.entry_points()
    return [ep for ep in eps.select(group="console_scripts")
            if (ep.dist is not None and ep.dist.name == "ligandparam")]


def test_console_scripts_are_declared():
    names = {ep.name for ep in _console_scripts()}
    assert names == {"lighfix", "lig-getparam", "smiles-to-pdb", "lig-to-sage"}, names


@pytest.mark.parametrize("ep", _console_scripts(), ids=lambda ep: ep.name)
def test_console_script_target_resolves(ep):
    """Each `[project.scripts]` target must actually exist and be callable.

    `lighfix` pointed at `cli_lighfix:lighfix` while the function was named `ligfix`,
    so the installed command raised ImportError for every release up to 1.0.0.
    """
    func = ep.load()
    assert callable(func)


def test_version_matches_distribution_metadata():
    """__version__ is derived from metadata, so it cannot drift from pyproject.toml."""
    assert ligandparam.__version__ == importlib.metadata.version("ligandparam")
