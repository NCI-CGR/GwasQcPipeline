from pathlib import Path
from shutil import copytree
from subprocess import run
from textwrap import dedent

import pytest

from cgr_gwas_qc.testing import chdir


@pytest.mark.workflow
def test_call_rate_filters(tmp_path: Path, bed_working_dir: Path):
    """Run call rate filter module and check outputs.

    Make sure outputs for CR1 and CR2 are created.
    """
    copytree(bed_working_dir, tmp_path, dirs_exist_ok=True)
    snake = tmp_path / "Snakefile"
    snake.write_text(
        dedent(
            """\
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("entry_points.smk")
            include: cfg.rules("snp_filters.smk")
            include: cfg.rules("sample_filters.smk")

            rule all:
                input:
                    expand("plink_filter_call_rate_2/snps.{ext}", ext=["bed", "bim", "fam"])
            """
        )
    )

    with chdir(tmp_path):
        run(["snakemake", "-j1", "--use-conda", "--nocolor"], check=True)

        assert Path("plink_filter_call_rate_1/snps.bed").exists()
        assert Path("plink_filter_call_rate_1/snps.bim").exists()
        assert Path("plink_filter_call_rate_1/snps.fam").exists()

        assert Path("plink_filter_call_rate_2/snps.bed").exists()
        assert Path("plink_filter_call_rate_2/snps.bim").exists()
        assert Path("plink_filter_call_rate_2/snps.fam").exists()
