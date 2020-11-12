from pathlib import Path
from shutil import copytree
from subprocess import run
from textwrap import dedent

import pytest

from cgr_gwas_qc.testing import chdir


@pytest.mark.workflow
def test_idat_intensity(tmp_path: Path, gtc_working_dir: Path):
    """Run idat_intensity tests."""
    copytree(gtc_working_dir, tmp_path, dirs_exist_ok=True)
    snake = tmp_path / "Snakefile"
    snake.write_text(
        dedent(
            """\
        from cgr_gwas_qc import load_config

        cfg = load_config()

        include: cfg.rules("common.smk")
        include: cfg.rules("entry_points.smk")
        include: cfg.rules("call_rate_filters.smk")
        include: cfg.rules("sample_filters.smk")

        rule all:
            input:
                "all_sample_idat_intensity/idat_intensity.csv",
        """
        )
    )

    with chdir(tmp_path):
        run(["snakemake", "-j1", "--use-conda", "--nocolor"], check=True)
        assert Path("all_sample_idat_intensity/idat_intensity.csv").exists()
