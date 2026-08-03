"""Tests for ligandparam.utils -- pure text/file manipulation, no external tools."""

import pytest

from ligandparam.utils import find_word_and_get_line, modify_gaussian_com


TWO_LINE = "%NPROC=16\n%MEM=32GB\n"


def _write(tmp_path, name, content):
    path = tmp_path / name
    path.write_text(content)
    return path


class TestModifyGaussianCom:
    """Gaussian documents Link 0 directives one per line, and G16 Rev. C.01 rejects
    the comma-separated single-line form. Both forms must be accepted on input and
    the two-line form must always be written."""

    def test_rewrites_legacy_one_line_header(self, tmp_path):
        com = _write(tmp_path, "legacy.com", "%NPROC=4, %MEM=8GB\n%chk=foo.chk\n#P HF/6-31G*\n")
        assert modify_gaussian_com(com, 16, 32) is True
        assert com.read_text() == TWO_LINE + "%chk=foo.chk\n#P HF/6-31G*\n"

    def test_rewrites_two_line_header(self, tmp_path):
        com = _write(tmp_path, "new.com", "%NPROC=4\n%MEM=8GB\n%chk=foo.chk\n")
        assert modify_gaussian_com(com, 16, 32) is True
        assert com.read_text() == TWO_LINE + "%chk=foo.chk\n"

    def test_accepts_crlf_header(self, tmp_path):
        """A CRLF-terminated header is matched and split across two lines.

        Only the header itself is rewritten, so the file keeps its original CRLF
        terminators elsewhere. Gaussian accepts the mixed endings.
        """
        com = tmp_path / "crlf.com"
        com.write_bytes(b"%NPROC=4\r\n%MEM=8GB\r\n%chk=foo.chk\r\n")
        assert modify_gaussian_com(com, 16, 32) is True
        assert com.read_bytes().startswith(b"%NPROC=16\n%MEM=32GB")

    def test_rewrites_every_link1_block(self, tmp_path):
        com = _write(
            tmp_path, "multi.com",
            "%NPROC=4, %MEM=8GB\n#P x\n--Link1--\n%NPROC=4, %MEM=8GB\n#P y\n")
        assert modify_gaussian_com(com, 16, 32) is True
        text = com.read_text()
        assert text.count("%NPROC=16\n%MEM=32GB") == 2
        assert "--Link1--" in text

    def test_returns_false_and_leaves_file_alone_without_a_header(self, tmp_path):
        original = "#P HF/6-31G*\n\ntitle\n\n0 1\n"
        com = _write(tmp_path, "noheader.com", original)
        assert modify_gaussian_com(com, 16, 32) is False
        assert com.read_text() == original


class TestFindWordAndGetLine:
    def test_finds_a_match_and_the_following_line(self, tmp_path):
        f = _write(tmp_path, "a.log", "one\nWarning! bad thing\ncontext line\nthree\n")
        assert find_word_and_get_line(f, "Warning!") == ["Warning! bad thing\ncontext line"]

    def test_match_on_last_line_without_trailing_newline(self, tmp_path):
        """The offending text must survive; this used to come back as ['']."""
        f = _write(tmp_path, "b.log", "line one\nWarning! at the very end")
        assert find_word_and_get_line(f, "Warning!") == ["Warning! at the very end"]

    def test_match_on_last_line_with_trailing_newline(self, tmp_path):
        f = _write(tmp_path, "c.log", "line one\nATTN needs attention\n")
        assert find_word_and_get_line(f, "ATTN") == ["ATTN needs attention"]

    def test_match_on_first_line(self, tmp_path):
        f = _write(tmp_path, "d.log", "ATTN right away\nsecond\n")
        assert find_word_and_get_line(f, "ATTN") == ["ATTN right away\nsecond"]

    def test_no_match_returns_empty_list(self, tmp_path):
        f = _write(tmp_path, "e.log", "nothing to see here\n")
        assert find_word_and_get_line(f, "ATTN") == []

    def test_empty_file_returns_empty_list(self, tmp_path):
        """mmap() rejects a zero-byte file; parmchk hits this when a run produces nothing."""
        f = _write(tmp_path, "empty.log", "")
        assert find_word_and_get_line(f, "ATTN") == []

    def test_multiple_matches(self, tmp_path):
        f = _write(tmp_path, "f.log", "ATTN one\nx\nATTN two\ny\n")
        assert len(find_word_and_get_line(f, "ATTN")) == 2
