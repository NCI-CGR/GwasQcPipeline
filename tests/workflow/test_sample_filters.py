from pathlib import Path
from shutil import copytree
from subprocess import run
from textwrap import dedent

import pytest

from cgr_gwas_qc.testing import chdir


@pytest.mark.workflow
def test_filter_contamination(tmp_path: Path, gtc_working_dir: Path):
    """Test sample contamination filter."""
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
                "sample_filters/not_contaminated/samples.bed",
                "sample_filters/not_contaminated/samples.bim",
                "sample_filters/not_contaminated/samples.fam",
        """
        )
    )

    with chdir(tmp_path):
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
