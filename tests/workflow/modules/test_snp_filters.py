from pathlib import Path
from subprocess import run
from textwrap import dedent

import pytest

from cgr_gwas_qc.testing import chdir, make_snakefile, make_test_config


################################################################################
# MAF Filtering
################################################################################
@pytest.mark.real_data
@pytest.fixture
def maf_filter_setup(conda_envs, real_data_cache, tmp_path) -> Path:
    conda_envs.copy_all_envs(tmp_path)

    real_data_cache.copy("original_data/manifest_full.csv", tmp_path / "manifest_full.csv")
    real_data_cache.copy(
        "production_outputs/plink_filter_call_rate_2", tmp_path / "plink_filter_call_rate_2"
    )

    make_test_config(tmp_path, data_type="real", sample_sheet="manifest_full.csv")

    make_snakefile(
        tmp_path,
        """
        from cgr_gwas_qc import load_config


        cfg = load_config()

        include: cfg.rules("common.smk")
        include: cfg.rules("snp_filters.smk")

        rule all:
            input:
                expand("snp_filters/maf/samples.{ext}", ext=["bed", "bim", "fam", "nosex"])
        """,
    )

    return tmp_path


@pytest.mark.regression
@pytest.mark.real_data
@pytest.mark.workflow
def test_maf_filter(maf_filter_setup):
    # GIVEN
    # - Snakefile, manifest.csv, conda-envs
    # - plink_filter_call_rate_2/samples.bed
    # - plink_filter_call_rate_2/samples.bim
    # - plink_filter_call_rate_2/samples.fam

    # WHEN run snakemake
    with chdir(maf_filter_setup):
        run(
            ["snakemake", "-j1", "--use-conda", "--nocolor", "--conda-frontend", "mamba"],
            check=True,
        )

        # THEN
        assert (maf_filter_setup / "snp_filters/maf/samples.bed").exists()
        assert (maf_filter_setup / "snp_filters/maf/samples.bim").exists()
        assert (maf_filter_setup / "snp_filters/maf/samples.fam").exists()
        assert (maf_filter_setup / "snp_filters/maf/samples.nosex").exists()


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
