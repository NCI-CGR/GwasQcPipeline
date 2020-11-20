from pathlib import Path
from shutil import copytree
from subprocess import run
from textwrap import dedent

import pytest

from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.testing.conda import CondaEnv


@pytest.mark.xfail
@pytest.mark.workflow
def test_snp_filters(tmp_path: Path, bed_working_dir: Path):
    """Test sample contamination filter."""
    CondaEnv().copy_env("plink2", tmp_path)
    copytree(bed_working_dir, tmp_path, dirs_exist_ok=True)
    snake = tmp_path / "Snakefile"
    snake.write_text(
        dedent(
            """\
        from cgr_gwas_qc import load_config

        cfg = load_config()

        include: cfg.rules("common.smk")
        include: cfg.rules("entry_points.smk")
        include: cfg.rules("call_rate_filters.smk")
        include: cfg.rules("snp_filters.smk")

        rule all:
            input:
                "snp_filters/ld_prune/samples.bed",
                "snp_filters/ld_prune/samples.bim",
                "snp_filters/ld_prune/samples.fam",
        """
        )
    )

    with chdir(tmp_path):
        run(
            ["snakemake", "-j1", "--use-conda", "--nocolor", "--conda-frontend", "mamba"],
            check=True,
        )

        assert Path("snp_filters/maf/samples.bed").exists()
        assert Path("snp_filters/maf/samples.bim").exists()
        assert Path("snp_filters/maf/samples.fam").exists()

        assert Path("snp_filters/ld_prune/ldPruneList.prune.in").exists()
        assert Path("snp_filters/ld_prune/ldPruneList.prune.out").exists()

        assert Path("snp_filters/ld_prune/samples.bed").exists()
        assert Path("snp_filters/ld_prune/samples.bim").exists()
        assert Path("snp_filters/ld_prune/samples.fam").exists()
