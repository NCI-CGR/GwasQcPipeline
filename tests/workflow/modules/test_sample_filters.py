from pathlib import Path
from subprocess import run
from textwrap import dedent

import pytest

from cgr_gwas_qc.testing import chdir


@pytest.mark.workflow
def test_filter_contamination(small_gtc_working_dir: Path, conda_envs):
    """Test sample contamination filter."""

    conda_envs.copy_all_envs(small_gtc_working_dir)
    snake = small_gtc_working_dir / "Snakefile"
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
                "sample_filters/not_contaminated/samples.bed",
                "sample_filters/not_contaminated/samples.bim",
                "sample_filters/not_contaminated/samples.fam",
        """
        )
    )

    with chdir(small_gtc_working_dir):
        run(
            ["snakemake", "-j1", "--use-conda", "--nocolor", "--conda-frontend", "mamba"],
            check=True,
        )
        assert Path("sample_filters/agg_median_idat_intensity.csv").exists()
        assert Path("sample_filters/small_manifest.AF.abf.txt").exists()
        assert Path("sample_filters/agg_contamination_test.csv").exists()
        assert Path("sample_filters/contaminated_samples.txt").exists()
        assert Path("sample_filters/not_contaminated/samples.bed").exists()
        assert Path("sample_filters/not_contaminated/samples.bim").exists()
        assert Path("sample_filters/not_contaminated/samples.fam").exists()


@pytest.mark.workflow
def test_ibd(tmp_path: Path, bed_working_dir: Path):
    """Test sample IBD."""
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
        include: cfg.rules("sample_filters.smk")

        rule all:
            input:
                "sample_filters/ibd/samples.genome",
        """
        )
    )

    with chdir(tmp_path):
        run(
            ["snakemake", "-j1", "--use-conda", "--nocolor", "--conda-frontend", "mamba"],
            check=True,
        )
        assert Path("sample_filters/ibd/samples.genome").exists()
