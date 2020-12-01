import os
from pathlib import Path
from textwrap import dedent

from cgr_gwas_qc import load_config
from cgr_gwas_qc.testing import (
    chdir,
    file_hashes_equal,
    file_rows_almost_equal,
    make_snakefile,
    make_test_config,
)


def test_chdir(tmp_path):
    """Test that chdir actually changes directory."""
    expected_in_chdir = Path(tmp_path).absolute().as_posix()

    before_chdir = Path(os.curdir).absolute().as_posix()
    with chdir(tmp_path):
        in_chdir = Path(os.curdir).absolute().as_posix()
    after_chdir = Path(os.curdir).absolute().as_posix()

    assert before_chdir == after_chdir
    assert expected_in_chdir == in_chdir
    assert in_chdir != before_chdir


def test_file_hashes_equal(tmp_path):
    """Identical files have the same hash."""
    file1 = tmp_path / "file1.txt"
    file2 = tmp_path / "file2.txt"

    file1.write_text("test")
    file2.write_text("test")

    assert file_hashes_equal(file1, file2)


def test_file_hashes_not_equal(tmp_path):
    """Different files do not have the same hash."""
    file1 = tmp_path / "file1.txt"
    file2 = tmp_path / "file2.txt"

    file1.write_text("test")
    file2.write_text("foo bar")

    assert file_hashes_equal(file1, file2) is False


def test_files_almost_equal(tmp_path):
    """Files with small float differences are the same."""
    file1 = tmp_path / "file1.txt"
    file2 = tmp_path / "file2.txt"

    file1.write_text("foo\t0.01\nbar\t0.05")
    file2.write_text("foo\t0.010001\nbar\t0.05")

    assert file_rows_almost_equal(file1, file2, fuzzy_col=1)


def test_files_almost_not_equal(tmp_path):
    """Files with big float differences are not the same."""
    file1 = tmp_path / "file1.txt"
    file2 = tmp_path / "file2.txt"

    file1.write_text("foo\t0.01\nbar\t0.05")
    file2.write_text("foo\t0.015\nbar\t0.05")

    assert file_rows_almost_equal(file1, file2, fuzzy_col=1) is False


def test_make_test_config(tmp_path):
    make_test_config(tmp_path)

    with chdir(tmp_path):
        cfg = load_config(validate=False)

    assert cfg.config.sample_sheet == "sample_sheet.csv"
    assert cfg.config.user_files.gtc_pattern is not None
    assert cfg.config.user_files.idat_pattern is not None


def test_make_snakefile(tmp_path):
    """Snakefile is created."""
    data = """
    from cgr_gwas_qc import load_config


    cfg = load_config()

    rule all:
        input: "test.txt"
    """

    make_snakefile(tmp_path, data)

    assert (tmp_path / "Snakefile").exists()
    assert dedent(data) == (tmp_path / "Snakefile").read_text()
