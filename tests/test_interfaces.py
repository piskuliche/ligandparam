"""Tests for the external-program interface layer.

These exercise the argument/redirection handling using /bin/cat and /bin/echo rather
than antechamber or Gaussian, so they need no AmberTools install.
"""

import os

import pytest

from ligandparam.interfaces import Gaussian, SimpleInterface


class Cat(SimpleInterface):
    def __init__(self, cwd, **kwargs):
        self.cwd = cwd
        self.logger = kwargs.get("logger") or __import__("logging").getLogger("test")
        self.set_method("cat")


def test_redirection_actually_pipes(tmp_path):
    """Passing a list with shell=True runs only argv[0] and drops every argument,
    so piped calls used to execute the bare binary and report success."""
    (tmp_path / "in.txt").write_text("hello world\n")
    Cat(tmp_path).call(inp_pipe="in.txt", out_pipe="out.txt")
    assert (tmp_path / "out.txt").read_text() == "hello world\n"


def test_redirection_resolves_relative_to_cwd(tmp_path):
    sub = tmp_path / "work"
    sub.mkdir()
    (sub / "in.txt").write_text("data\n")
    Cat(sub).call(inp_pipe="in.txt", out_pipe="out.txt")
    assert (sub / "out.txt").exists()
    assert not (tmp_path / "out.txt").exists()


def test_absolute_redirection_paths_are_honoured(tmp_path):
    src = tmp_path / "abs_in.txt"
    src.write_text("abs\n")
    dst = tmp_path / "abs_out.txt"
    Cat(tmp_path).call(inp_pipe=str(src), out_pipe=str(dst))
    assert dst.read_text() == "abs\n"


def test_missing_input_file_raises_filenotfound(tmp_path):
    """More precise than the previous RuntimeError-with-shell-text; RuntimeError is
    reserved for a non-zero exit code from the program itself."""
    with pytest.raises(FileNotFoundError):
        Cat(tmp_path).call(inp_pipe="does_not_exist.txt", out_pipe="out.txt")


def test_nonzero_exit_code_raises_runtimeerror(tmp_path):
    failing = Cat(tmp_path)
    failing.set_method("false")
    with pytest.raises(RuntimeError):
        failing.call()


def test_dry_run_does_not_execute(tmp_path):
    (tmp_path / "in.txt").write_text("hello\n")
    Cat(tmp_path).call(inp_pipe="in.txt", out_pipe="out.txt", dry_run=True)
    assert not (tmp_path / "out.txt").exists()


def test_nproc_does_not_leak_into_the_parent_environment(tmp_path):
    """`env = os.environ` is a reference: setting OMP_NUM_THREADS on it pinned every
    later stage in the process to the first stage's thread count."""
    before = os.environ.get("OMP_NUM_THREADS")
    iface = Cat(tmp_path)
    iface.nproc = 1
    (tmp_path / "in.txt").write_text("x\n")
    iface.call(inp_pipe="in.txt", out_pipe="out.txt")
    assert os.environ.get("OMP_NUM_THREADS") == before


class TestGaussianDefaults:
    def test_binary_defaults_to_g16(self, tmp_path):
        """kwargs.get(opt, "") made the documented default unreachable, and an empty
        binary produced a submit script that was a bare redirect."""
        assert Gaussian(cwd=tmp_path).method == "g16"

    def test_explicit_binary_is_kept(self, tmp_path):
        assert Gaussian(cwd=tmp_path, gaussian_binary="g09").method == "g09"

    def test_empty_binary_falls_back(self, tmp_path):
        assert Gaussian(cwd=tmp_path, gaussian_binary="").method == "g16"

    def test_missing_cwd_is_rejected(self):
        with pytest.raises(ValueError):
            Gaussian()

    def test_environment_is_not_mutated(self, tmp_path):
        before = os.environ.get("GAUSS_SCRDIR")
        g = Gaussian(cwd=tmp_path, gaussian_scratch="/scratch/xyz")
        env = g.set_environment()
        assert env["GAUSS_SCRDIR"] == "/scratch/xyz"
        assert os.environ.get("GAUSS_SCRDIR") == before
