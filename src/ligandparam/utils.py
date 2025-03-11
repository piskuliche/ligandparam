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

def modify_gaussian_com(filepath: Path, nproc:int, mem:int):
    config_line_regex = re.compile(b"%NPROC=\d+, %MEM=\d+MB")
    nproc_bytes = str(nproc).encode()
    mem_bytes = str(mem).encode()
    new_line_prefix = b"%NPROC="
    new_line_sep = b", %MEM="
    new_line_suffix = b"MB"
    new_line = new_line_prefix + nproc_bytes + new_line_sep + mem_bytes + new_line_suffix

    with open(filepath, 'r+b') as f:
        mm = mmap.mmap(f.fileno(), 0)
        try:
            output_parts = []
            last_match_end = 0
            for match in config_line_regex.finditer(mm):
                start_index = match.start()
                end_index = match.end()
                output_parts.append(mm[last_match_end:start_index])
                output_parts.append(new_line)
                last_match_end = end_index
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
    return False