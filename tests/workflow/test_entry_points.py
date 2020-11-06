from pathlib import Path
from shutil import copytree
from subprocess import run

import pytest

from cgr_gwas_qc.testing import chdir

snakefile = """\
from cgr_gwas_qc import load_config

cfg = load_config()

include: cfg.rules("entry_points.smk")

rule all:
    input:
        expand("plink_start/samples.{ext}", ext=["bed", "bim", "fam"])
"""


def run_and_check_for_plink_start_samples(tmp_path):
    with chdir(tmp_path):
        run(["snakemake", "-j1", "--use-conda", "--nocolor"], check=True)
        assert Path("plink_start/samples.bed").exists()
        assert Path("plink_start/samples.bim").exists()
        assert Path("plink_start/samples.fam").exists()


@pytest.mark.workflow
def test_gtc_entry(tmp_path: Path, gtc_working_dir: Path):
    copytree(gtc_working_dir, tmp_path, dirs_exist_ok=True)
    snake = tmp_path / "Snakefile"
    snake.write_text(snakefile)
    run_and_check_for_plink_start_samples(tmp_path)


@pytest.mark.workflow
def test_ped_entry(tmp_path: Path, ped_working_dir: Path):
    copytree(ped_working_dir, tmp_path, dirs_exist_ok=True)
    snake = tmp_path / "Snakefile"
    snake.write_text(snakefile)
    run_and_check_for_plink_start_samples(tmp_path)


@pytest.mark.workflow
def test_bed_entry(tmp_path: Path, bed_working_dir: Path):
    copytree(bed_working_dir, tmp_path, dirs_exist_ok=True)
    snake = tmp_path / "Snakefile"
    snake.write_text(snakefile)
    run_and_check_for_plink_start_samples(tmp_path)

    # Check files are symbolic links
    assert Path("plink_start/samples.bed").is_symlink()
    assert Path("plink_start/samples.bim").is_symlink()
    assert Path("plink_start/samples.fam").is_symlink()
