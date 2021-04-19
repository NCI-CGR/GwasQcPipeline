import shutil

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from cgr_gwas_qc.testing import file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.workflow
@pytest.mark.regression
@pytest.mark.real_data
def test_sample_qc_table(sample_qc_df):
    # GIVEN: The sample qc report
    # THEN: This should be identical to the production outputs except for:
    mapper = {
        "Call_Rate_1": {"name": "Call_Rate_1", "dtype": "float"},
        "Call_Rate_2": {"name": "Call_Rate_2", "dtype": "float"},
        "Call_Rate_Initial": {"name": "Call_Rate_Initial", "dtype": "float"},
        "ChrX_Inbreed_estimate": {"name": "X_inbreeding_coefficient", "dtype": "float"},
        "Contamination_Rate": {"name": "Contamination_Rate", "dtype": "float"},
        "Count_of_QC_Issue": {"name": "Count_of_QC_Issue", "dtype": "UInt8"},
        "Expected Replicate Discordance": {"name": "is_replicate_discordant", "dtype": "boolean"},
        "IdatIntensity": {"name": "IdatIntensity", "dtype": "float"},
        "Identifiler_Needed": {"name": "identifiler_needed", "dtype": "boolean"},
        "Low Call Rate": {"name": "is_call_rate_filtered", "dtype": "boolean"},
        "Sample_ID": {"name": "Sample_ID", "dtype": "string"},
        "Sex Discordant": {"name": "is_sex_discordant", "dtype": "boolean"},
        "Unexpected Replicate": {"name": "is_unexpected_replicate", "dtype": "boolean"},
    }

    exp_df = (
        pd.read_csv(RealData() / "production_outputs/all_sample_qc.csv")
        .rename({k: v["name"] for k, v in mapper.items()}, axis=1)
        .reindex([v["name"] for v in mapper.values()], axis=1)
        .astype({v["name"]: v["dtype"] for v in mapper.values()})
        .set_index("Sample_ID")
        .sort_index()
    )

    obs_df = (
        sample_qc_df.reindex([v["name"] for v in mapper.values()], axis=1)
        .set_index("Sample_ID")
        .sort_index()
        .fillna({"is_call_rate_filtered": True})  # Legacy will have True instead of pd.NA
    )

    assert_frame_equal(exp_df, obs_df)


@pytest.mark.workflow
@pytest.mark.real_data
def test_sample_qc_stats(tmp_path, sample_qc_csv):
    """Test the createion of summary_stats.txt

    I am not able to perform regression testing because there are some issues
    with the original script. First, they have the contamination threshold of
    0.1 hardcoded so the summary stats don't match the sample_qc_table if
    this threshold is changed by the user. Second, the fail call rate
    includes NAs, where in the sample_qc_table these are converted to True.
    Third, because I switched from R to python, there are some small
    naming/formatting differences that I felt were not necessary to repeat.
    """
    # GIVEN: real data sample sheet, config, and sample_qc_summary table.
    (
        RealData(tmp_path)
        .make_cgr_sample_sheet()
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("sample_level_summary.smk")

            rule all:
                input:
                    "sample_level/sample_qc_summary_stats.txt"
            """
        )
    )
    (tmp_path / "sample_level").mkdir()
    shutil.copy(sample_qc_csv, tmp_path / "sample_level/sample_qc.csv")

    # WHEN: run snakemake to create sample_qc_summary_stats.txt
    run_snakemake(tmp_path)

    # THEN: The file should end with my sample counts
    assert (tmp_path / "sample_level/sample_qc_summary_stats.txt").read_text().endswith("203  17\n")


@pytest.mark.workflow
@pytest.mark.regression
@pytest.mark.real_data
def test_qc_failures(tmp_path, sample_qc_csv):
    # GIVEN: real data sample sheet, config, and sample_qc_summary table.
    data_cache = (
        RealData(tmp_path)
        .make_cgr_sample_sheet()
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("sample_level_summary.smk")

            rule all:
                input:
                    "sample_level/qc_failures/low_call_rate.txt",
                    "sample_level/qc_failures/contaminated.txt",
                    "sample_level/qc_failures/sex_discordant.txt",
                    "sample_level/qc_failures/replicate_discordant.txt",
                    "sample_level/internal_controls.txt",
            """
        )
    )
    (tmp_path / "sample_level").mkdir()
    shutil.copy(sample_qc_csv, tmp_path / "sample_level/sample_qc.csv")

    # WHEN: run snakemake to create sample_qc_summary_stats.txt
    run_snakemake(tmp_path)

    # THEN: The files should match production outputs
    obs_ = tmp_path / "sample_level/qc_failures/low_call_rate.txt"
    exp_ = data_cache / "production_outputs/remove_qc_fail/LowCR.txt"
    file_hashes_equal(obs_, exp_)

    obs_ = tmp_path / "sample_level/qc_failures/contaminated.txt"
    exp_ = data_cache / "production_outputs/remove_qc_fail/contaminated.txt"
    file_hashes_equal(obs_, exp_)

    obs_ = tmp_path / "sample_level/qc_failures/sex_discordant.txt"
    exp_ = data_cache / "production_outputs/remove_qc_fail/SexDiscordant.txt"
    file_hashes_equal(obs_, exp_)

    obs_ = tmp_path / "sample_level/qc_failures/replicate_discordant.txt"
    exp_ = data_cache / "production_outputs/remove_qc_fail/ExpectedRepDiscordant.txt"
    file_hashes_equal(obs_, exp_)

    obs_ = tmp_path / "sample_level/internal_controls.txt"
    exp_ = data_cache / "production_outputs/subject_level/internalControls.txt"
    file_hashes_equal(obs_, exp_)


@pytest.mark.workflow
@pytest.mark.regression
@pytest.mark.real_data
def test_remove_contaminated(tmp_path, conda_envs):
    # GIVEN: real data sample sheet, config
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .make_cgr_sample_sheet()
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples.bed",
            "sample_level/call_rate_2/samples.bed",
        )
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples.bim",
            "sample_level/call_rate_2/samples.bim",
        )
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples.fam",
            "sample_level/call_rate_2/samples.fam",
        )
        .copy(
            "production_outputs/remove_qc_fail/contaminated.txt",
            "sample_level/qc_failures/contaminated.txt",
        )
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("sample_level_summary.smk")

            rule all:
                input:
                    "sample_level/call_rate_2/samples_contaminated_removed.bed"
            """
        )
    )

    # WHEN:
    run_snakemake(tmp_path)

    # THEN:
    obs_ = tmp_path / "sample_level/call_rate_2/samples_contaminated_removed.bed"
    exp_ = data_cache / "production_outputs/remove_qc_fail/samples.bed"
    assert file_hashes_equal(obs_, exp_)

    obs_ = tmp_path / "sample_level/call_rate_2/samples_contaminated_removed.bim"
    exp_ = data_cache / "production_outputs/remove_qc_fail/samples.bim"
    assert file_hashes_equal(obs_, exp_)

    obs_ = tmp_path / "sample_level/call_rate_2/samples_contaminated_removed.fam"
    exp_ = data_cache / "production_outputs/remove_qc_fail/samples.fam"
    assert file_hashes_equal(obs_, exp_)


@pytest.mark.workflow
@pytest.mark.real_data
def test_remove_contaminated_add_contam_samples(tmp_path, conda_envs):
    """Run removal step forcing 3 contaminated samples.

    Unfortunately, none of my 206 test samples are actually flagged as
    contaminated. Here I am using the Sex Discordant samples instead to force
    the `remove_contaminated` rule to actually do something.
    """
    # GIVEN: real data sample sheet, config, 3 "contaminated" samples
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .make_cgr_sample_sheet()
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples.bed",
            "sample_level/call_rate_2/samples.bed",
        )
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples.bim",
            "sample_level/call_rate_2/samples.bim",
        )
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples.fam",
            "sample_level/call_rate_2/samples.fam",
        )
        .copy(
            "production_outputs/remove_qc_fail/SexDiscordant.txt",
            "sample_level/qc_failures/contaminated.txt",
        )  # I am calling the 3 sex discordant samples as contaminated
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("sample_level_summary.smk")

            rule all:
                input:
                    "sample_level/call_rate_2/samples_contaminated_removed.bed"
            """
        )
    )

    # WHEN:
    run_snakemake(tmp_path)

    # THEN: The output files should not match the legacy run
    obs_ = tmp_path / "sample_level/call_rate_2/samples_contaminated_removed.bed"
    exp_ = data_cache / "production_outputs/remove_qc_fail/samples.bed"
    assert not file_hashes_equal(obs_, exp_)

    obs_ = tmp_path / "sample_level/call_rate_2/samples_contaminated_removed.bim"
    exp_ = data_cache / "production_outputs/remove_qc_fail/samples.bim"
    assert not file_hashes_equal(obs_, exp_)

    obs_ = tmp_path / "sample_level/call_rate_2/samples_contaminated_removed.fam"
    exp_ = data_cache / "production_outputs/remove_qc_fail/samples.fam"
    assert not file_hashes_equal(obs_, exp_)

    # The plink log should say we remove 3 samples
    assert (
        "--remove: 203 people"
        in (tmp_path / "sample_level/call_rate_2/samples_contaminated_removed.log").read_text()
    )
