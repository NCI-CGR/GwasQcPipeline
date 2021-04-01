from pathlib import Path

import pandas as pd
import pytest

from cgr_gwas_qc.testing.data import RealData


@pytest.mark.real_data
def test_ArrayProcessing(real_config, sample_sheet_full):
    from cgr_gwas_qc.reporting.sample_qc import ArrayProcessing

    ap = ArrayProcessing.construct(real_config, sample_sheet_full)

    assert 220 == ap.num_samples_qc_processed
    assert 0 == ap.num_failed_array_processing


@pytest.mark.real_data
def test_CompletionRate(snp_qc_df, sample_qc_df):
    from cgr_gwas_qc.reporting.sample_qc import CompletionRate

    cr = CompletionRate.construct(snp_qc_df, sample_qc_df, Path("test.png"))

    assert pytest.approx(0.978, rel=0.001) == cr.mean_initial_call_rate

    assert 3738 == cr.num_snp_cr1_filtered  # NOTE: Legacy has 5218 b/c it includes NaNs
    assert 694860 == cr.num_snp_pass_cr1

    assert 4 == cr.num_sample_cr1_filtered
    assert 216 == cr.num_sample_pass_cr1

    assert 12666 == cr.num_snp_cr2_filtered
    assert 682194 == cr.num_snp_pass_cr2

    assert 10 == cr.num_sample_cr2_filtered
    assert 206 == cr.num_sample_pass_cr2

    assert 206 == cr.num_sample_pass_call_rate
    assert 682194 == cr.num_snp_pass_call_rate


@pytest.mark.real_data
def test_Contamination(sample_qc_df):
    from cgr_gwas_qc.reporting.sample_qc import Contamination

    cn = Contamination.construct(sample_qc_df)

    assert 0 == cn.num_pass_cr_and_contaminated
    assert 206 == cn.num_remaining


@pytest.mark.real_data
@pytest.fixture
def control_replicates_df(pytestconfig) -> pd.DataFrame:
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    return pd.read_csv(RealData() / "production_outputs/concordance/InternalQcKnown.csv").rename(
        {"Concordance": "concordance"}, axis=1
    )


@pytest.mark.real_data
def test_InternalControls(sample_qc_df, control_replicates_df):
    from cgr_gwas_qc.reporting.sample_qc import InternalControls

    ic = InternalControls.construct(sample_qc_df, control_replicates_df)

    assert 7 == ic.num_internal_controls
    assert 199 == ic.num_remaining
    assert pytest.approx(0.9999, rel=1e-4) == ic.min_concordance
    assert pytest.approx(0.9999, rel=1e-4) == ic.mean_concordance


@pytest.mark.real_data
@pytest.fixture
def study_replicates_df(pytestconfig) -> pd.DataFrame:
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    return pd.read_csv(RealData() / "production_outputs/concordance/StudySampleKnown.csv").rename(
        {"Concordance": "concordance"}, axis=1
    )


@pytest.mark.real_data
def test_ExpectedReplicates(sample_qc_df, study_replicates_df):
    from cgr_gwas_qc.reporting.sample_qc import ExpectedReplicates

    er = ExpectedReplicates.construct(sample_qc_df, study_replicates_df)

    assert 0 == er.num_low_concordance
    assert 15 == er.num_not_subject_representative
    assert 184 == er.num_remaining
    assert pytest.approx(0.9996, rel=1e-4) == er.min_concordance
    assert pytest.approx(0.9999, rel=1e-4) == er.mean_concordance


@pytest.mark.real_data
@pytest.fixture
def unexpected_replicates_df(pytestconfig) -> pd.DataFrame:
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    return pd.read_csv(RealData() / "production_outputs/concordance/UnknownReplicates.csv").rename(
        {"Concordance": "concordance"}, axis=1
    )


@pytest.mark.real_data
def test_UnExpectedReplicates(sample_qc_df, unexpected_replicates_df):
    from cgr_gwas_qc.reporting.sample_qc import UnExpectedReplicates

    unexpected = UnExpectedReplicates.construct(sample_qc_df, unexpected_replicates_df)

    assert 0 == unexpected.num_unexpected_replicates


@pytest.mark.real_data
def test_SampleSummary(sample_qc_df):
    from cgr_gwas_qc.reporting.sample_qc import SampleSummary

    sm = SampleSummary.construct(sample_qc_df)

    assert "67:146" == sm.case_control
    assert "127:93" == sm.male_female
    assert 197 == sm.num_starting_subjects
    assert 36 == sm.num_samples_excluded
    assert 184 == sm.num_subjects_remaining


@pytest.mark.real_data
def test_SexVerification(sample_qc_df):
    from cgr_gwas_qc.reporting.sample_qc import SexVerification

    sv = SexVerification.construct(sample_qc_df, Path("test.png"))

    assert 3 == sv.num_sex_discordant


@pytest.mark.real_data
def test_Ancestry(sample_qc_df):
    from cgr_gwas_qc.reporting.sample_qc import Ancestry

    ac = Ancestry.construct(sample_qc_df, Path("test.png"))

    expected = (
        "| Ancestry Group   |   Count |\n"
        "|:-----------------|--------:|\n"
        "| ADMIXED_AFR      |       1 |\n"
        "| AFR              |       8 |\n"
        "| AFR_EUR          |       7 |\n"
        "| ASN              |       5 |\n"
        "| EUR              |     163 |"
    )

    assert expected == ac.table
