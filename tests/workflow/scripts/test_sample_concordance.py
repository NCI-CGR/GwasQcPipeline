import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from cgr_gwas_qc.reporting.constants import REPORT_NAME_MAPPER
from cgr_gwas_qc.testing.data import RealData
from cgr_gwas_qc.workflow.scripts import concordance_table, sample_concordance


@pytest.mark.real_data
@pytest.fixture
def sample_concordance_table_df(sample_concordance_table_csv) -> pd.DataFrame:
    return concordance_table.read(sample_concordance_table_csv).rename(
        {"ID1": "Sample_ID1", "ID2": "Sample_ID2"}, axis=1
    )


@pytest.mark.real_data
@pytest.fixture
def known_replicates_df(real_cfg, sample_concordance_table_df):
    return sample_concordance._known_replicates_df(
        sample_concordance_table_df,
        sample_concordance._get_known_replicates(real_cfg.ss),
        real_cfg.ss.set_index("Sample_ID").Group_By_Subject_ID.rename("Subject_ID"),
    )


@pytest.mark.real_data
@pytest.mark.regression
def test_known_replicates_df(known_replicates_df):
    # GIVEN: Expected output from legacy workflow
    exp_df = pd.read_csv(RealData() / "production_outputs/concordance/KnownReplicates.csv")
    exp_tweak = exp_df.set_index(["Subject_ID", "Sample_ID1", "Sample_ID2"]).reindex(
        ["Concordance", "PI_HAT"], axis=1
    )

    # Observed output from dev workflow
    obs_df = known_replicates_df
    obs_tweak = (
        obs_df.set_index(["Subject_ID", "Sample_ID1", "Sample_ID2"])
        .rename(REPORT_NAME_MAPPER, axis=1)
        .reindex(["Concordance", "PI_HAT"], axis=1)
    )

    # THEN: Frames should be the same.
    assert_frame_equal(exp_tweak, obs_tweak, check_dtype=False, check_exact=False, check_like=True)


@pytest.mark.real_data
def test_known_replicates_df_force_no_replicates(real_cfg, sample_concordance_table_df):
    """Test what happens when there are no known replicates."""
    # GIVEN: That we have no known replicates
    obs_df = sample_concordance._known_replicates_df(
        sample_concordance_table_df,
        [],  # force no known replicates
        real_cfg.ss.set_index("Sample_ID").Group_By_Subject_ID.rename("Subject_ID"),
    )

    # THEN: The results table should be empty
    assert 0 == obs_df.shape[0]


@pytest.mark.real_data
@pytest.mark.regression
def test_unknown_replicates_df(real_cfg, sample_concordance_table_df):
    """The test data set has no unknown replicates in it."""
    # GIVEN: Expected output from legacy workflow
    exp_df = pd.read_csv(RealData() / "production_outputs/concordance/UnknownReplicates.csv")
    exp_tweak = exp_df.set_index(
        ["Subject_ID1", "Subject_ID2", "Sample_ID1", "Sample_ID2"]
    ).reindex(["Concordance", "PI_HAT"], axis=1)

    # Observed output from dev workflow
    obs_df = sample_concordance._unknown_replicates_df(
        sample_concordance_table_df,
        sample_concordance._get_known_replicates(real_cfg.ss),
        real_cfg.ss.set_index("Sample_ID").Group_By_Subject_ID.rename("Subject_ID"),
    )
    obs_tweak = (
        obs_df.set_index(["Subject_ID1", "Subject_ID2", "Sample_ID1", "Sample_ID2"])
        .rename(REPORT_NAME_MAPPER, axis=1)
        .reindex(["Concordance", "PI_HAT"], axis=1)
    )

    # THEN: Frames should be the same.
    assert_frame_equal(exp_tweak, obs_tweak, check_dtype=False, check_exact=False, check_like=True)


@pytest.mark.real_data
def test_unknown_replicates_df_force_replicates(real_cfg, sample_concordance_table_df):
    """Force unknown replicates.

    Since the legacy workflow had no unknown replicates, I want to force the
    table to be built. I am doing this by just saying there are no "known"
    replicates.
    """
    # GIVEN: That I say there are no known replicates
    obs_df = sample_concordance._unknown_replicates_df(
        sample_concordance_table_df,
        [],  # Force no known replicates
        real_cfg.ss.set_index("Sample_ID").Group_By_Subject_ID.rename("Subject_ID"),
    )

    # THEN: all concordant samples should be added to this table.
    assert sample_concordance_table_df.is_ge_concordance.sum() == obs_df.shape[0]


@pytest.mark.real_data
def test_split_into_qc(real_cfg, known_replicates_df):
    # GIVEN: Expected output from legacy workflow
    exp_df = pd.read_csv(RealData() / "production_outputs/concordance/InternalQcKnown.csv")
    exp_tweak = exp_df.set_index(["Subject_ID", "Sample_ID1", "Sample_ID2"]).reindex(
        ["Concordance", "PI_HAT"], axis=1
    )

    # WHEN: I split the known replicates into qc fraction
    qc, _ = sample_concordance._split_into_qc_and_study_samples(real_cfg.ss, known_replicates_df)

    qc_tweak = (
        qc.set_index(["Subject_ID", "Sample_ID1", "Sample_ID2"])
        .rename(REPORT_NAME_MAPPER, axis=1)
        .reindex(["Concordance", "PI_HAT"], axis=1)
    )

    # THEN: Frames should be the same.
    assert_frame_equal(exp_tweak, qc_tweak, check_dtype=False, check_exact=False, check_like=True)


@pytest.mark.real_data
def test_split_into_study(real_cfg, known_replicates_df):
    # GIVEN: Expected output from legacy workflow
    exp_df = pd.read_csv(RealData() / "production_outputs/concordance/StudySampleKnown.csv")
    exp_tweak = exp_df.set_index(["Subject_ID", "Sample_ID1", "Sample_ID2"]).reindex(
        ["Concordance", "PI_HAT"], axis=1
    )

    # WHEN: I split the known replicates into study fraction
    _, study = sample_concordance._split_into_qc_and_study_samples(real_cfg.ss, known_replicates_df)

    study_tweak = (
        study.set_index(["Subject_ID", "Sample_ID1", "Sample_ID2"])
        .rename(REPORT_NAME_MAPPER, axis=1)
        .reindex(["Concordance", "PI_HAT"], axis=1)
    )

    # THEN: Frames should be the same.
    assert_frame_equal(
        exp_tweak, study_tweak, check_dtype=False, check_exact=False, check_like=True
    )
