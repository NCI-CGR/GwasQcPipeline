import pandas as pd
import pytest
from pandas.testing import assert_frame_equal, assert_series_equal

from cgr_gwas_qc.testing import run_snakemake
from cgr_gwas_qc.testing.data import RealData


################################################################################
# Legacy Regression Tests
################################################################################
@pytest.mark.real_data
@pytest.mark.regression
def test_legacy_agg_median_idat_intensity(real_data_cache):
    legacy = (
        pd.read_csv(real_data_cache / "legacy_outputs/all_sample_idat_intensity/idat_intensity.csv")
        .rename(
            {
                "SampId": "Sample_ID",
                "ChipId": "Chip_ID",
                "MedianIntensity": "median_intensity",
            },
            axis=1,
        )
        .set_index("Sample_ID")
    )

    dev = pd.read_csv(
        real_data_cache / "dev_outputs/sample_level/contamination/median_idat_intensity.csv"
    ).set_index("Sample_ID")

    assert_frame_equal(legacy, dev, check_like=True)


@pytest.mark.real_data
@pytest.mark.regression
def test_legacy_agg_contamination(real_data_cache):
    legacy = (
        pd.read_csv(real_data_cache / "legacy_outputs/all_contam/contam.csv")
        .rename({"ID": "Sample_ID"}, axis=1)
        .set_index("Sample_ID")
    )

    dev = pd.read_csv(
        real_data_cache / "dev_outputs/sample_level/contamination/verifyIDintensity.csv"
    ).set_index("Sample_ID")

    mask = legacy["%Mix"].notna() & dev["%Mix"].notna()

    assert_series_equal(legacy[mask]["%Mix"], dev[mask]["%Mix"], atol=1e-3)
    assert_series_equal(legacy[mask]["LLK"], dev[mask]["LLK"], atol=1000)
    assert_series_equal(legacy[mask]["LLK0"], dev[mask]["LLK0"], atol=500)


@pytest.mark.real_data
@pytest.mark.regression
def test_legacy_abf(real_data_cache):
    legacy_abf = real_data_cache / "legacy_outputs/GSAMD-24v1-0_20011747_A1.AF.abf.txt"
    dev_abf = real_data_cache / "dev_outputs/sample_level/GSAMD-24v1-0_20011747_A1.AF.abf.txt"

    legacy = pd.read_csv(legacy_abf, sep="\t")
    dev = pd.read_csv(dev_abf, sep="\t").fillna(0)  # legacy uses 0.0 instead of NA.

    # There were some logic changes so these file are moderately different.
    prop_very_different = (abs(legacy.ABF - dev.ABF) > 1e-6).mean()
    assert 0.01 > prop_very_different  # less than 1% are more than 0.000001 different


################################################################################
# Workflow Tests
################################################################################
@pytest.mark.slow
@pytest.mark.real_data
@pytest.fixture(scope="module")
def contamination(pytestconfig, tmp_path_factory, conda_envs):
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    tmp_path = tmp_path_factory.mktemp("contamination")
    conda_envs.copy_env("illuminaio", tmp_path)
    conda_envs.copy_env("verifyidintensity", tmp_path)

    data_cache = (
        RealData(tmp_path, full_sample_sheet=False)
        .add_user_files(entry_point="gtc", copy=False)
        .make_config(num_snps=700078)
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            module contam:
                snakefile: cfg.subworkflow("contamination")
                config: {}

            use rule * from contam
            """
        )
    )

    # Copy ABF file b/c it takes a long time to create
    data_cache.copy(
        "dev_outputs/sample_level/GSAMD-24v1-0_20011747_A1.AF.abf.txt",
        "sample_level/GSAMD-24v1-0_20011747_A1.AF.abf.txt",
    )

    run_snakemake(tmp_path, keep_temp=True)

    return tmp_path


@pytest.mark.slow
@pytest.mark.workflow
@pytest.mark.real_data
def test_agg_median_idat_intensity(real_data_cache, contamination):
    dev = pd.read_csv(
        real_data_cache / "dev_outputs/sample_level/contamination/median_idat_intensity.csv"
    ).set_index("Sample_ID")

    snake = (
        pd.read_csv(contamination / "sample_level/contamination/median_idat_intensity.csv")
        .set_index("Sample_ID")
        .astype({"median_intensity": float})
    )

    assert_frame_equal(dev.reindex(snake.index), snake)


@pytest.mark.slow
@pytest.mark.workflow
@pytest.mark.real_data
def test_agg_contamination_test(real_data_cache, contamination):
    dev = pd.read_csv(
        real_data_cache / "dev_outputs/sample_level/contamination/verifyIDintensity.csv"
    ).set_index("Sample_ID")

    snake = pd.read_csv(
        contamination / "sample_level/contamination/verifyIDintensity.csv"
    ).set_index("Sample_ID")

    assert_frame_equal(dev.reindex(snake.index), snake)
