from hashlib import sha256
from math import isclose
from pathlib import Path

from cgr_gwas_qc.typing import PathLike


################################################################################
# General
################################################################################
def file_hashes_equal(file1: PathLike, file2: PathLike) -> bool:
    """Comapres two files using sha256 hashes.

    Takes two files and calculates their hashes returning True if they are
    the same and False otherwise.
    """
    return sha256(Path(file1).read_bytes()).digest() == sha256(Path(file2).read_bytes()).digest()


def sorted_file_equal(file1: PathLike, file2: PathLike) -> bool:
    """Comapres two sorted files.

    Reads in each file, splits by line, and then sorts.
    """
    pth1, pth2 = Path(file1), Path(file2)
    return sorted(pth1.read_text().splitlines()) == sorted(pth2.read_text().splitlines())


def file_rows_almost_equal(
    file1: PathLike,
    file2: PathLike,
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


################################################################################
# PLINK File Fromats
################################################################################
def assert_plink_bed_equal(file1: PathLike, file2: PathLike):
    raise NotImplementedError


def assert_plink_bim_equal(file1: PathLike, file2: PathLike):
    """Compare two PLINK BIM files ignoring allele order.

    Raises:
        AssertionError: if all rows are not equal after ignoring allele order.
    """
    from cgr_gwas_qc.parsers.bim import BimFile

    fh1, fh2 = BimFile(file1), BimFile(file2)
    rows_are_equal = [record1 == record2 for record1, record2 in zip(fh1, fh2)]
    assert all(rows_are_equal)


def assert_plink_fam_equal(file1: PathLike, file2: PathLike):
    """Compare two PLINK FAM files.

    Raises:
        AssertionError: if files are not identical.
    """
    assert file_hashes_equal(file1, file2)
