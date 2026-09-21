"""Low-level file and stream helpers used across ligandparam."""

import ctypes
import io
import mmap
import os
import re
import sys
import tempfile
from contextlib import contextmanager

libc = None
c_stderr = None
if sys.platform == "darwin":
    libc = ctypes.CDLL(None)
    c_stderr = ctypes.c_void_p.in_dll(libc, "__stderrp")
elif sys.platform.startswith("linux"):
    libc = ctypes.CDLL(None)
    c_stderr = ctypes.c_void_p.in_dll(libc, "stderr")
else:
    # Windows and other platforms: C-level stderr redirect is unavailable.
    pass

from pathlib import Path
from typing import Optional,  Union


def find_word_and_get_line(filepath: Union[Path, str], word: str):
    """Return lines containing ``word`` using memory-mapped search.

    Parameters
    ----------
    filepath : path-like
        File to search.
    word : str
        Substring to find.

    Returns
    -------
    list of str
        Matching lines (stripped).
    """
    word_b = word.encode()  # Encode the word to bytes for searching in mmap
    lines_found = []

    with open(filepath, mode='rb') as file:
        with mmap.mmap(file.fileno(), length=0, access=mmap.ACCESS_READ) as mm:
            for match in re.finditer(word_b, mm):  # Use re.finditer to find all occurrences
                start = match.start()
                end = match.end()

                # Find the start of the line (go back until newline or start of file)
                line_start = mm.rfind(b'\n', 0, start) + 1  # +1 to move past the newline
                if line_start == -1:
                    line_start = 0  # Handle case where match is on the first line

                # Find the end of the line and then get the next line as well.
                line_end = mm.find(b'\n', end)
                line_end = mm.find(b'\n', line_end + 1)
                if line_end == -1:
                    line_end = mm.size()  # Handle case where match is on the last line

                # Extract and decode the line
                line = mm[line_start:line_end].decode('utf-8')  # Adjust decoding if needed
                lines_found.append(line.strip())

    return lines_found


@contextmanager
def stderr_redirector(stream):
    """Temporarily redirect C-level stderr into ``stream``.

    Parameters
    ----------
    stream : io.IOBase
        Destination stream for captured stderr bytes.

    Yields
    ------
    None
    """
    # The original fd stderr points to. Usually 1 on POSIX systems.
    original_stderr_fd = sys.stderr.fileno()

    def _redirect_stderr(to_fd):
        """Redirect stderr to the given file descriptor."""
        # Flush the C-level buffer stderr
        if libc is not None and c_stderr is not None:
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
        # Yield to caller, then redirect stderr back to the saved fd
        yield
        _redirect_stderr(saved_stderr_fd)
        # Copy contents of temporary file to the given stream
        tfile.flush()
        tfile.seek(0, io.SEEK_SET)
        stream.write(tfile.read())
    finally:
        tfile.close()
        os.close(saved_stderr_fd)
