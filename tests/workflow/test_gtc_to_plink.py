from pathlib import Path
from shutil import copytree
from subprocess import run
from textwrap import dedent

import pytest

from cgr_gwas_qc.testing import chdir


@pytest.mark.workflow
def test_gtc_to_plink(tmp_path: Path, working_dir: Path):
    copytree(working_dir, tmp_path, dirs_exist_ok=True)
    snake = tmp_path / "Snakefile"
    snake.write_text(
        dedent(
            """\
    from cgr_gwas_qc import load_config

    cfg = load_config()

    include: cfg.rules("entry_points.smk")

    rule all:
        input:
            expand("plink_start/samples.{ext}", ext=["bed", "bim", "fam"])
    """
        )
    )

    with chdir(tmp_path):
        run(["snakemake", "-j1", "--use-conda", "--nocolor"], check=True)
        assert Path("plink_start/samples.bed").exists()
        assert Path("plink_start/samples.bim").exists()
        assert Path("plink_start/samples.fam").exists()
