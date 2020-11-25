"""Tests for the entry_points module.

Each entry point is expected to create the following files:

- plink_start/samples.bed
- plink_start/samples.bim
- plink_start/samples.ped

Entry point decisions are made by the presence of different user provided
files in the config. If these options are None (or not in the file) then a
different entry point is used.

For testing entry points we need to make sure that the entry point files are
present in the working directory and that they are referenced in the config.
"""
from pathlib import Path
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
    """Helper function to make sure expected outputs are present."""
    with chdir(tmp_path):
        run(["snakemake", "-j1", "--use-conda", "--nocolor"], check=True)
        assert Path("plink_start/samples.bed").exists()
        assert Path("plink_start/samples.bim").exists()
        assert Path("plink_start/samples.fam").exists()


@pytest.mark.workflow
def test_gtc_entry(small_gtc_working_dir: Path, conda_envs):
    """GTC entry point.

    This is the primary entry point and requires sample level GTC files.
    """
    conda_envs.copy_all_envs(small_gtc_working_dir)
    snake = small_gtc_working_dir / "Snakefile"
    snake.write_text(snakefile)
    run_and_check_for_plink_start_samples(small_gtc_working_dir)


@pytest.mark.workflow
def test_ped_entry(small_ped_working_dir: Path, conda_envs):
    """Aggregated PED entry point.

    This entry point expects an aggregated (multi-sample) PED/MAP file.
    """
    conda_envs.copy_all_envs(small_ped_working_dir)
    snake = small_ped_working_dir / "Snakefile"
    snake.write_text(snakefile)
    run_and_check_for_plink_start_samples(small_ped_working_dir)


@pytest.mark.workflow
def test_bed_entry(small_bed_working_dir: Path):
    """Aggregated BED entry point.

    This entry point expects an aggregated (multi-sample) BED/BIM/FAM file.
    Since this is the expected output of the entry point module, these files
    are just symbolically linked to the expected location
    ``plink_start/samples.{bed,bim,fam}``.
    """
    snake = small_bed_working_dir / "Snakefile"
    snake.write_text(snakefile)
    run_and_check_for_plink_start_samples(small_bed_working_dir)

    with chdir(small_bed_working_dir):
        # Check files are symbolic links
        assert Path("plink_start/samples.bed").is_symlink()
        assert Path("plink_start/samples.bim").is_symlink()
        assert Path("plink_start/samples.fam").is_symlink()
