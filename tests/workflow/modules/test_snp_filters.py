from pathlib import Path
from subprocess import run
from textwrap import dedent

import pytest

from cgr_gwas_qc.testing import chdir


@pytest.mark.xfail
@pytest.mark.workflow
def test_snp_filters(small_bed_working_dir: Path, conda_envs):
    """Test sample contamination filter."""
    conda_envs.copy_all_envs(small_bed_working_dir)
    snake = small_bed_working_dir / "Snakefile"
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

    with chdir(small_bed_working_dir):
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
