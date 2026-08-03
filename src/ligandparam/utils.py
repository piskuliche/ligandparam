import ctypes
import io
import mmap
import os
import re
import sys
import tempfile
import warnings
from contextlib import contextmanager

libc = ctypes.CDLL(None)
if sys.platform == "darwin":
    c_stderr = ctypes.c_void_p.in_dll(libc, "__stderrp")
elif sys.platform.startswith("linux"):
    c_stderr = ctypes.c_void_p.in_dll(libc, "stderr")
else:
    # Leaving c_stderr undefined turned this into a NameError deep inside
    # stderr_redirector(). None makes libc.fflush() flush all open streams instead.
    warnings.warn(f"Unsupported platform {sys.platform}: cannot look up the C stderr "
                  f"symbol, so stderr redirection will flush all streams instead.")
    c_stderr = None


from pathlib import Path
from typing import Optional,  Union


def find_word_and_get_line(filepath: Union[Path, str], word: str):
    """
    Finds a word in a file using memory mapping and returns the full lines
    containing the word.

    Parameters
    ----------
    filepath : Union[Path, str]
        The path to the file to search.
    word : str
        The word to search for in the file.

    Returns
    -------
    list of str
        A list of lines containing the word.
    """
    word_b = word.encode()  # Encode the word to bytes for searching in mmap
    lines_found = []

    if os.path.getsize(filepath) == 0:
        # mmap() rejects empty files, and an empty file trivially contains nothing.
        return lines_found

    with open(filepath, mode='rb') as file:
        with mmap.mmap(file.fileno(), length=0, access=mmap.ACCESS_READ) as mm:
            for match in re.finditer(word_b, mm):  # Use re.finditer to find all occurrences
                start = match.start()
                end = match.end()

                # Find the start of the line. rfind returns -1 when the match is on the
                # first line, and -1 + 1 == 0 is already the answer we want.
                line_start = mm.rfind(b'\n', 0, start) + 1

                # Take the matching line plus the following one, stopping at EOF. Each
                # find must be guarded: searching from -1 + 1 == 0 would otherwise wrap
                # back to the top of the file and produce line_end < line_start, which
                # yields an empty string.
                line_end = mm.find(b'\n', end)
                if line_end == -1:
                    line_end = mm.size()
                else:
                    next_end = mm.find(b'\n', line_end + 1)
                    line_end = mm.size() if next_end == -1 else next_end

                # Extract and decode the line
                line = mm[line_start:line_end].decode('utf-8')  # Adjust decoding if needed
                lines_found.append(line.strip())

    return lines_found


def modify_gaussian_com(filepath: Path, nproc: int, mem: int):
    """
    Modifies a Gaussian input file to update the number of processors and memory allocation.

    Parameters
    ----------
    filepath : Path
        The path to the Gaussian input file.
    nproc : int
        The number of processors to set in the file.
    mem : int
        The amount of memory (in GB) to set in the file.

    Returns
    -------
    bool
        True if at least one Link 0 header was rewritten, False if the file contained
        no recognisable header (in which case it is left untouched).
    """
    # Match either the legacy one-line form or the one-directive-per-line form that
    # Gaussian documents (and that G16 Rev. C.01 requires), but always write the latter.
    config_line_regex = re.compile(rb"%NPROC=\d+(?:, |\r?\n)%MEM=\d+GB")
    nproc_bytes = str(nproc).encode()
    mem_bytes = str(mem).encode()
    new_line_prefix = b"%NPROC="
    new_line_sep = b"\n%MEM="
    new_line_suffix = b"GB"
    new_line = new_line_prefix + nproc_bytes + new_line_sep + mem_bytes + new_line_suffix

    with open(filepath, 'r+b') as f:
        mm = mmap.mmap(f.fileno(), 0)
        try:
            output_parts = []
            last_match_end = 0
            nmatches = 0
            for match in config_line_regex.finditer(mm):
                start_index = match.start()
                end_index = match.end()
                output_parts.append(mm[last_match_end:start_index])
                output_parts.append(new_line)
                last_match_end = end_index
                nmatches += 1
            if nmatches == 0:
                return False
            output_parts.append(mm[last_match_end:])
            new_content = b"".join(output_parts)

            mm.close()
            f.seek(0)
            f.truncate(0)
            f.write(new_content)
            return True
        finally:
            if 'mm' in locals() and not mm.closed:
                mm.close()


@contextmanager
def stderr_redirector(stream):
    """
    Redirects stderr to a given stream within a context.

    Parameters
    ----------
    stream : io.IOBase
        The stream to which stderr will be redirected.

    Yields
    ------
    None
        Allows the caller to execute code with stderr redirected.
    """
    # The original fd stderr points to. Usually 1 on POSIX systems.
    original_stderr_fd = sys.stderr.fileno()

    def _redirect_stderr(to_fd):
        """Redirect stderr to the given file descriptor."""
        # Flush the C-level buffer stderr
        libc.fflush(c_stderr)
        # Flush and close sys.stderr - also closes the file descriptor (fd)
        sys.stderr.close()
        # Make original_stderr_fd point to the same file as to_fd
        os.dup2(to_fd, original_stderr_fd)
        # Create a new sys.stderr that points to the redirected fd
        sys.stderr = io.TextIOWrapper(os.fdopen(original_stderr_fd, 'wb'))

    # Save a copy of the original stderr fd in saved_stderr_fd
    saved_stderr_fd = os.dup(original_stderr_fd)
    # Create a temporary file and redirect stderr to it
    tfile = tempfile.TemporaryFile(mode='w+b')
    try:
        _redirect_stderr(tfile.fileno())
        yield
    finally:
        # The restore has to happen even when the caller raises. Leaving it after the
        # yield meant an exception inside the block left fd 2 pointing at a temp file
        # that `finally` then closed, so the traceback for the original error vanished.
        try:
            _redirect_stderr(saved_stderr_fd)
            # Copy contents of temporary file to the given stream
            tfile.flush()
            tfile.seek(0, io.SEEK_SET)
            stream.write(tfile.read())
        finally:
            tfile.close()
            os.close(saved_stderr_fd)
