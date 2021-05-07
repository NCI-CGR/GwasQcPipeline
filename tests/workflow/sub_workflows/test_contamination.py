import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_series_equal

from cgr_gwas_qc.parsers.illumina import AdpcReader
from cgr_gwas_qc.testing import file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.real_data
@pytest.mark.workflow
@pytest.mark.slow
@pytest.fixture(scope="module")
def contam_output(tmp_path_factory, conda_envs):
    tmp_path = tmp_path_factory.mktemp("contamination_workflow")
    conda_envs.copy_env("illuminaio", tmp_path)
    conda_envs.copy_env("verifyidintensity", tmp_path)
    (
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

            use rule * from contam as contam_*
            """
        )
    )

    # WHEN: I run snakemake to calculate median IDAT intensity and aggregate
    # these results.
    run_snakemake(tmp_path, keep_temp=True)

    return tmp_path


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_agg_median_idat_intensity(contam_output):
    # GIVEN: aggregated median idat intensity from the contamination workflow
    filename = contam_output / "sample_level/median_idat_intensity.csv"

    # THEN: The values in the observed summary table match the values found in
    # the legacy files.
    obs_table = (
        pd.read_csv(filename)  # Sample_ID Chip_ID median_intensity
        .set_index("Chip_ID")
        .median_intensity
    )
    for exp_file in (RealData() / "production_outputs/idat_intensity").glob("*.txt"):
        exp_value = int(exp_file.read_text())
        exp_chip_id = exp_file.stem.split(".")[0]
        assert obs_table[exp_chip_id] == exp_value


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_per_sample_gtc_to_adpc(contam_output):
    # GIVEN: contamination subworkflow output
    ss = pd.read_csv(contam_output / "cgr_sample_sheet.csv")

    # THEN: The binary adpc.bin files should be the same as production. There
    # may be slight difference between floating point number
    for r in ss.itertuples(index=False):
        obs_adpc = contam_output / f"sample_level/per_sample_adpc//{r.Sample_ID}.adpc.bin"
        exp_adpc = RealData() / f"production_outputs/contam/{r.Sample_ID}.adpc.bin"
        with AdpcReader(obs_adpc) as obs, AdpcReader(exp_adpc) as expect:
            for exp_row, obs_row in zip(expect, obs):
                # each int should be the exact same
                assert exp_row.x_raw == obs_row.x_raw
                assert exp_row.y_raw == obs_row.y_raw
                assert exp_row.genotype == obs_row.genotype

                if np.isnan(obs_row.x_norm) | np.isnan(obs_row.y_norm):
                    # Note: the legacy workflow outputs (i.e., exp_row) will
                    # have 0.0 instead of NA for x_norm and y_norm. This
                    # difference is okay, see Issue #31 for details.

                    # If normalized values are NA then the following should be true
                    assert (obs_row.x_raw, obs_row.y_raw) == (0, 0)
                    assert obs_row.genotype_score == 0.0
                    assert obs_row.genotype == 3
                else:
                    # floats should be really close
                    assert pytest.approx(exp_row.x_norm, abs=1e-6) == obs_row.x_norm
                    assert pytest.approx(exp_row.y_norm, abs=1e-6) == obs_row.y_norm
                    assert pytest.approx(exp_row.genotype_score, abs=1e-6) == obs_row.genotype_score


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.slow
def test_per_sample_verifyIDintensity_contamination(tmp_path, conda_envs):
    # GIVEN: Real data with the per sample ADPC.bin files, the 1KG B allele
    # frequencies, and verifyIDintensity conda environment.
    conda_envs.copy_env("verifyidintensity", tmp_path)
    data_cache = (
        RealData(tmp_path, full_sample_sheet=False)
        .add_user_files(entry_point="gtc", copy=False)
        .copy("production_outputs/contam", "sample_level/per_sample_adpc")
        .copy(
            "production_outputs/GSAMD-24v1-0_20011747_A1.AF.abf.txt",
            "sample_level/GSAMD-24v1-0_20011747_A1.AF.abf.txt",
        )
        .make_config(num_snps=700078)
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config
            cfg = load_config()

            include: cfg.subworkflow("contamination")

            rule all:
                input:
                    cfg.expand("sample_level/per_sample_contamination_test/{Sample_ID}.contam.out")
            """
        )
    )

    # WHEN: I run snakemake to create per sample contamination estimates
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: The observed file should be identical to the expected file for each
    # sample.
    for r in data_cache.ss.data.itertuples(index=False):
        obs_file = tmp_path / f"sample_level/per_sample_contamination_test/{r.Sample_ID}.contam.out"
        exp_file = (
            data_cache / f"production_outputs/one_samp_b_1000g_contam/{r.Sample_ID}.contam.out"
        )
        assert file_hashes_equal(obs_file, exp_file)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_contamination_test_with_missing_abf_values(contam_output):
    """Does verifyIDintensity run with missing abf values.

    In issue #30 we found that the legacy pipeline was setting missing values
    to 0.0 when converting BPM to abf. After looking at verifyIDintensity
    code we found that these values should be `"NA"` to be excluded from the
    contamination calculation. Here I am testing that everything runs and
    generates similar outputs.
    """
    # GIVEN: Real data using GTC entry point and the per sample adpc files.
    # GIVEN: contamination subworkflow output
    ss = pd.read_csv(contam_output / "cgr_sample_sheet.csv")

    # THEN: The observed contamination file is similar to the expected file for each sample.
    for r in ss.itertuples(index=False):
        obs_file = (
            contam_output / f"sample_level/per_sample_contamination_test/{r.Sample_ID}.contam.out"
        )
        exp_file = (
            RealData() / f"production_outputs/one_samp_b_1000g_contam/{r.Sample_ID}.contam.out"
        )

        obs = obs_file.read_text().strip().split("\n")[-1].split()
        exp_ = exp_file.read_text().strip().split("\n")[-1].split()
        assert pytest.approx(float(exp_[1]), abs=0.01) == float(obs[1])  # within 1% of each other
        assert pytest.approx(float(exp_[2]), rel=0.1) == float(obs[2])  # within 1000
        assert pytest.approx(float(exp_[3]), rel=0.1) == float(obs[3])  # within 1000


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.slow
def test_contamination_test_with_missing_adpc_values(tmp_path, conda_envs):
    """Does verifyIDintensity work when there are missing data.

    In Issue #31 I determined that there is a small difference when building
    the adpc files. The current scripts outputs normalized values as `nan`
    when raw values are 0. The legacy outputs them as 0. I want to make sure
    that verifyIDintensity does not change regardless of `nan` or 0.
    """
    # GIVEN: Real data with the 1KG B allele frequencies and verifyIDintensity
    # conda environment.
    conda_envs.copy_env("verifyidintensity", tmp_path)
    data_store = (
        RealData(tmp_path, full_sample_sheet=False)
        .add_user_files(entry_point="gtc", copy=False)
        .copy(
            "production_outputs/GSAMD-24v1-0_20011747_A1.AF.abf.txt",
            "sample_level/GSAMD-24v1-0_20011747_A1.AF.abf.txt",
        )
        .make_config(num_snps=700078)
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("common")
            include: cfg.subworkflow("contamination")

            rule all:
                input:
                    cfg.expand("sample_level/per_sample_contamination_test/{Sample_ID}.contam.out")
            """
        )
    )

    # WHEN: I run snakemake to create per sample contamination estimates
    # (adpc.bin) and contamination estimates (contam.out)
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: I created the adbc.bin files
    for r in data_store.ss.data.itertuples(index=False):
        assert (tmp_path / f"sample_level/per_sample_adpc/{r.Sample_ID}.adpc.bin").exists()

    # The observed contamination file is identical to the expected file for each sample.
    for r in data_store.ss.data.itertuples(index=False):
        obs_file = tmp_path / f"sample_level/per_sample_contamination_test/{r.Sample_ID}.contam.out"
        exp_file = (
            data_store / f"production_outputs/one_samp_b_1000g_contam/{r.Sample_ID}.contam.out"
        )
        assert file_hashes_equal(obs_file, exp_file)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_agg_contamination_test(contam_output):

    # THEN: The observed and subset of the expected have almost identical
    # values.
    obs_df = (
        pd.read_csv(contam_output / "sample_level/contamination.csv")
        .set_index("Sample_ID")
        .drop("is_contaminated", axis=1)
        .dropna(how="any")
    )

    exp_df = (
        pd.read_csv(RealData() / "production_outputs/all_contam/contam.csv")
        .rename({"ID": "Sample_ID"}, axis=1)
        .set_index("Sample_ID")
        .reindex(obs_df.index)
    )

    assert_series_equal(exp_df["%Mix"], obs_df["%Mix"], atol=1e-4)
    assert_series_equal(exp_df["LLK"], obs_df["LLK"], rtol=0.05)
    assert_series_equal(exp_df["LLK0"], obs_df["LLK0"], rtol=0.05)
