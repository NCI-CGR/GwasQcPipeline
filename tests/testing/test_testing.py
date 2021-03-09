import os
from pathlib import Path
from textwrap import dedent

from cgr_gwas_qc.testing import chdir, file_hashes_equal, file_rows_almost_equal, make_snakefile


def test_chdir(tmp_path):
    """Test that chdir actually changes directory."""
    # GIVEN: our current working directory (start_dir) and a temporary directory (tmp_dir)
    start_dir = Path(os.curdir).absolute().as_posix()
    tmp_dir = Path(tmp_path).absolute().as_posix()

    # WHEN:
    with chdir(tmp_dir):
        # we enter the context manager our working directory changes to tmp_dir
        in_chdir = Path(os.curdir).absolute().as_posix()
    # we leave the context manager our working directory goes back to start_dir
    after_chdir = Path(os.curdir).absolute().as_posix()

    # THEN:
    # our start_dir and tmp_dir are different
    assert start_dir != tmp_dir
    # after leaving the context manager start_dir and working dir are the same
    assert start_dir == after_chdir
    # in the context manager tmp_dir and the working directory are the same
    assert tmp_dir == in_chdir


def test_file_hashes_equal(tmp_path):
    """Identical files have the same hash."""
    # GIVEN: two identical files
    file1 = tmp_path / "file1.txt"
    file1.write_text("test")

    file2 = tmp_path / "file2.txt"
    file2.write_text("test")

    # WHEN-THEN: we compare their hashes they are the same
    assert file_hashes_equal(file1, file2)


def test_file_hashes_not_equal(tmp_path):
    """Different files do not have the same hash."""
    # GIVEN: two different files
    file1 = tmp_path / "file1.txt"
    file1.write_text("test")

    file2 = tmp_path / "file2.txt"
    file2.write_text("foo bar")

    # WHEN-THEN: we compare their hashes they are different
    assert file_hashes_equal(file1, file2) is False


def test_files_almost_equal(tmp_path):
    """Files with small float differences are the same."""
    # GIVEN: two files with a column of very similar floating point numbers
    file1 = tmp_path / "file1.txt"
    file1.write_text("foo\t0.01\nbar\t0.05")

    file2 = tmp_path / "file2.txt"
    file2.write_text("foo\t0.010001\nbar\t0.05")

    # WHEN-THEN: we compare them they are the same
    assert file_rows_almost_equal(file1, file2, fuzzy_col=1)


def test_files_almost_not_equal(tmp_path):
    """Files with big float differences are not the same."""
    # GIVEN: two files with a column of slightly similar floating point numbers
    file1 = tmp_path / "file1.txt"
    file1.write_text("foo\t0.01\nbar\t0.05")

    file2 = tmp_path / "file2.txt"
    file2.write_text("foo\t0.015\nbar\t0.05")

    # WHEN-THEN: we compare them they are the different
    assert file_rows_almost_equal(file1, file2, fuzzy_col=1) is False


def test_make_snakefile(tmp_path):
    """Snakefile is created."""
    # GIVEN: a string representation of a Snakefile and tmp_path
    data = """
    from cgr_gwas_qc import load_config


    cfg = load_config()

    rule all:
        input: "test.txt"
    """

    # WHEN: write the string to a Snakefile in tmp_path
    make_snakefile(tmp_path, data)

    # THEN:
    # the `Snakefile` exists
    assert (tmp_path / "Snakefile").exists()
    # the contents are the same as our string version
    assert dedent(data) == (tmp_path / "Snakefile").read_text()
