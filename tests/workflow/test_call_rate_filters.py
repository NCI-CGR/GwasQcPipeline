from pathlib import Path
from shutil import copytree
from subprocess import run
from textwrap import dedent

import pytest

from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.testing.conda import CondaEnv


@pytest.mark.workflow
def test_call_rate_filters(tmp_path: Path, bed_working_dir: Path):
    """Run call rate filter module and check outputs.

    Make sure outputs for CR1 and CR2 are created.
    """
    CondaEnv().copy_env("plink2", tmp_path)
    copytree(bed_working_dir, tmp_path, dirs_exist_ok=True)
    snake = tmp_path / "Snakefile"
    snake.write_text(
        dedent(
            """\
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("entry_points.smk")
            include: cfg.rules("call_rate_filters.smk")

            rule all:
                input:
                    expand("plink_filter_call_rate_2/samples.{ext}", ext=["bed", "bim", "fam"])
            """
        )
    )

    with chdir(tmp_path):
        run(["snakemake", "-j1", "--use-conda", "--nocolor"], check=True)

        assert Path("plink_filter_call_rate_1/snps.bed").exists()
        assert Path("plink_filter_call_rate_1/snps.bim").exists()
        assert Path("plink_filter_call_rate_1/snps.fam").exists()
        assert Path("plink_filter_call_rate_1/snps.log").exists()

        assert Path("plink_filter_call_rate_1/samples.bed").exists()
        assert Path("plink_filter_call_rate_1/samples.bim").exists()
        assert Path("plink_filter_call_rate_1/samples.fam").exists()
        assert Path("plink_filter_call_rate_1/samples.log").exists()

        assert Path("plink_filter_call_rate_2/snps.bed").exists()
        assert Path("plink_filter_call_rate_2/snps.bim").exists()
        assert Path("plink_filter_call_rate_2/snps.fam").exists()
        assert Path("plink_filter_call_rate_2/snps.log").exists()

        assert Path("plink_filter_call_rate_2/samples.bed").exists()
        assert Path("plink_filter_call_rate_2/samples.bim").exists()
        assert Path("plink_filter_call_rate_2/samples.fam").exists()
        assert Path("plink_filter_call_rate_2/samples.log").exists()


@pytest.mark.parametrize(
    "prefix", ["plink_start", "plink_filter_call_rate_1", "plink_filter_call_rate_2"]
)
@pytest.mark.workflow
def test_call_rate_statistics(tmp_path: Path, bed_working_dir: Path, prefix: str):
    """Calculate call rate statistics."""
    CondaEnv().copy_env("plink2", tmp_path)
    copytree(bed_working_dir, tmp_path, dirs_exist_ok=True)
    snake = tmp_path / "Snakefile"
    snake.write_text(
        dedent(
            f"""\
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("entry_points.smk")
            include: cfg.rules("call_rate_filters.smk")

            rule all:
                input:
                    expand("{prefix}/samples.{{ext}}", ext=["imiss", "lmiss"])
            """
        )
    )

    with chdir(tmp_path):
        run(["snakemake", "-j1", "--use-conda", "--nocolor"], check=True)

        assert Path(f"{prefix}/samples.imiss").exists()
        assert Path(f"{prefix}/samples.lmiss").exists()
