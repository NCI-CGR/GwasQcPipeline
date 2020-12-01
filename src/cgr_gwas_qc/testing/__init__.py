import contextlib
import os
from hashlib import sha256
from math import isclose
from pathlib import Path
from textwrap import dedent
from typing import Union


@contextlib.contextmanager
def chdir(dirname: Path):
    curdir = Path().cwd()
    try:
        os.chdir(dirname)
        yield
    finally:
        os.chdir(curdir)


def file_hashes_equal(file1: Union[str, Path], file2: Union[str, Path]) -> bool:
    """Comapres two files using sha256 hashes.

    Takes two files and calculates their hashes returning True if they are
    the same and False otherwise.
    """
    return sha256(Path(file1).read_bytes()).digest() == sha256(Path(file2).read_bytes()).digest()


def file_rows_almost_equal(
    file1: Union[str, Path],
    file2: Union[str, Path],
    fuzzy_col: int,
    sep: str = "\t",
    header: bool = False,
) -> bool:
    """Compares two files row by row and makes sure they are almost equal"""
    results = []
    for row_idx, (r1, r2) in enumerate(
        zip(Path(file1).read_text().splitlines(), Path(file2).read_text().splitlines())
    ):
        if header and row_idx == 0:
            results.append(r1 == r2)
            continue

        for idx, (c1, c2) in enumerate(zip(r1.split(sep), r2.split(sep))):
            if idx != fuzzy_col:
                results.append(c1 == c2)
            else:
                results.append(isclose(float(c1), float(c2), rel_tol=1e-4))

    return all(results)


def make_snakefile(working_dir: Path, contents: str):
    snakefile = working_dir / "Snakefile"
    snakefile.write_text(dedent(contents))
