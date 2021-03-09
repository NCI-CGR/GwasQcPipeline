import shutil
from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from cgr_gwas_qc.testing import file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.workflow
@pytest.mark.real_data
@pytest.fixture(scope="session")
def sample_qc_report(tmp_path_factory) -> Path:
    # GIVEN: All sample QC results (including contamination). NOTE: I am using
    # SNPweights ancestry b/c I don't have production GRAF output
    tmp_path = tmp_path_factory.mktemp("sample_level_summary")
    (
        RealData(tmp_path)
        .copy_sample_sheet()
        .add_user_files(entry_point="gtc", copy=False)
        .copy("production_outputs/plink_start/samples_start.imiss", "sample_level/samples.imiss")
        .copy(
            "production_outputs/plink_filter_call_rate_1/samples_filter1.imiss",
            "sample_level/call_rate_1/samples.imiss",
        )
        .copy(
            "production_outputs/plink_filter_call_rate_1/samples_filter1.sexcheck",
            "sample_level/call_rate_1/samples.sexcheck",
        )
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss",
            "sample_level/call_rate_2/samples.imiss",
        )
        .copy(
            "production_outputs/snpweights/samples.snpweights.csv",
            "sample_level/ancestry/graf_ancestry.txt",
        )
        .copy(
            "production_outputs/concordance/KnownReplicates.csv",
            "sample_level/concordance/KnownReplicates.csv",
        )
        .copy(
            "production_outputs/concordance/UnknownReplicates.csv",
            "sample_level/concordance/UnknownReplicates.csv",
        )
        .copy(
            "production_outputs/all_contam/contam.csv",
            "sample_level/contamination/verifyIDintensity_contamination.csv",
        )
        .copy(
            "production_outputs/all_sample_idat_intensity/idat_intensity.csv",
            "sample_level/median_idat_intensity.csv",
        )
        .make_config(workflow_params={"subject_id_to_use": "PI_Subject_ID"})
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("sample_level_summary.smk")

            rule all:
                input:
                    "sample_level/qc_summary.csv"
            """
        )
    )

    # WHEN: I run snakemake to generate the sample_qc_report
    run_snakemake(tmp_path)

    return tmp_path / "sample_level/qc_summary.csv"


@pytest.mark.workflow
@pytest.mark.regression
@pytest.mark.real_data
def test_sample_qc_report(sample_qc_report):
    # GIVEN: The sample qc report
    # THEN: This should be identical to the production outputs except for:
    exclude_cols = [
        "IdatsInProjectDir",  # This column does not match b/c I did not have all the Idat files
        "Identifiler_Reason",  # New column
        "Internal_Control",  # New column
        "Group_By_Subject_ID",  # New column
        "Subject_Representative",  # New column
        "Subject_Dropped_From_Study",  # New column
        "Case/Control_Status",  # New column
    ]

    obs_ = (
        pd.read_csv(sample_qc_report)
        .drop(exclude_cols, axis=1, errors="ignore")
        .sort_values("Sample_ID")
        .reset_index(drop=True)
    )
    exp_ = (
        pd.read_csv(RealData() / "production_outputs/all_sample_qc.csv")
        .drop(exclude_cols, axis=1, errors="ignore")
        .sort_values("Sample_ID")
        .reset_index(drop=True)
    )

    assert_frame_equal(obs_, exp_)


@pytest.mark.workflow
@pytest.mark.real_data
def test_sample_qc_stats(tmp_path, sample_qc_report):
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
        .copy_sample_sheet()
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("sample_level_summary.smk")

            rule all:
                input:
                    "sample_level/qc_summary_stats.txt"
            """
        )
    )
    (tmp_path / "sample_level").mkdir()
    shutil.copy(sample_qc_report, tmp_path / "sample_level/qc_summary.csv")

    # WHEN: run snakemake to create qc_summary_stats.txt
    run_snakemake(tmp_path)

    # THEN: The file should end with my sample counts
    assert (tmp_path / "sample_level/qc_summary_stats.txt").read_text().endswith("203  17\n")


@pytest.mark.workflow
@pytest.mark.regression
@pytest.mark.real_data
def test_qc_failures(tmp_path, sample_qc_report):
    # GIVEN: real data sample sheet, config, and sample_qc_summary table.
    data_cache = (
        RealData(tmp_path)
        .copy_sample_sheet()
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
    shutil.copy(sample_qc_report, tmp_path / "sample_level/qc_summary.csv")

    # WHEN: run snakemake to create qc_summary_stats.txt
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
        .copy_sample_sheet()
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
        .copy_sample_sheet()
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
