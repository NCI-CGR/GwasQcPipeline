"""Tests for the call_rate_filters module.

This module uses plink to successively filter samples that fall below two
different call rate thresholds. Samples and SNPs are filtered using the same
``plink`` command, but I don't know if the filtering is always done in the
same order.

Outputs from this module includes:

- Call Rate 1
    - plink_filter_call_rate_1/samples.bed
    - plink_filter_call_rate_1/samples.bim
    - plink_filter_call_rate_1/samples.fam

- Call Rate 2
    - plink_filter_call_rate_2/samples.bed
    - plink_filter_call_rate_2/samples.bim
    - plink_filter_call_rate_2/samples.fam
"""
from pathlib import Path
from shutil import copytree
from subprocess import run

import pytest

from cgr_gwas_qc.testing import chdir

snakefile = """\
from cgr_gwas_qc import load_config

cfg = load_config()

include: cfg.rules("entry_points.smk")
include: cfg.rules("call_rate_filters.smk")

rule all:
    input:
        expand("plink_filter_call_rate_2/samples.{ext}", ext=["bed", "bim", "fam"])
"""


@pytest.mark.workflow
def test_call_rate_filters(tmp_path: Path, bed_working_dir: Path):
    """Run call rate filter module and check outputs.

    Make sure outputs for CR1 and CR2 are created.
    """
    copytree(bed_working_dir, tmp_path, dirs_exist_ok=True)
    snake = tmp_path / "Snakefile"
    snake.write_text(snakefile)

    with chdir(tmp_path):
        run(["snakemake", "-j1", "--use-conda", "--nocolor"], check=True)

        assert Path("plink_filter_call_rate_1/samples.bed").exists()
        assert Path("plink_filter_call_rate_1/samples.bim").exists()
        assert Path("plink_filter_call_rate_1/samples.fam").exists()

        assert Path("plink_filter_call_rate_2/samples.bed").exists()
        assert Path("plink_filter_call_rate_2/samples.bim").exists()
        assert Path("plink_filter_call_rate_2/samples.fam").exists()
