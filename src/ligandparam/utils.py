import mmap
import re

from pathlib import Path
from typing import Union


def find_word_and_get_line(filepath: Union[Path, str], word: str):
    """
    Finds a word in a file using memory mapping and returns the full lines
    containing the word.
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
                line_end = mm.find(b'\n', line_end+1)
                if line_end == -1:
                    line_end = mm.size()  # Handle case where match is on the last line

                # Extract and decode the line
                line = mm[line_start:line_end].decode('utf-8')  # Adjust decoding if needed
                lines_found.append(line.strip())

    return lines_found

