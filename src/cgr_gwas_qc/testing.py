from hashlib import sha256
from pathlib import Path
from typing import Union


def file_hashes_equal(file1: Union[str, Path], file2: Union[str, Path]) -> bool:
    """Comapres two files using sha256 hashes.

    Takes two files and calculates their hashes returning True if they are
    the same and False otherwise.
    """
    return sha256(Path(file1).read_bytes()).digest() == sha256(Path(file2).read_bytes()).digest()
