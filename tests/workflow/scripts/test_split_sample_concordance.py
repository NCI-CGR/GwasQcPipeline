import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from cgr_gwas_qc.workflow.scripts import split_sample_concordance


@pytest.mark.real_data
def test_known_full(split_sample_concordance_tables, real_data_cache):
    exp_df = pd.read_csv(
        real_data_cache / "production_outputs/concordance/KnownReplicates.csv"
    ).set_index(["Sample_ID1", "Sample_ID2"])

    obs_df = (
        split_sample_concordance.read_known_sample_concordance(
            split_sample_concordance_tables / "KnownReplicates.csv"
        )
        .set_index(["Sample_ID1", "Sample_ID2"])
        .rename(
            {
                "PLINK_PI_HAT": "PI_HAT",
                "PLINK_concordance": "Concordance",
                "Subject_ID1": "Subject_ID",
            },
            axis=1,
        )
        .reindex(exp_df.columns, axis=1)
    )

    assert_frame_equal(exp_df, obs_df, check_like=True, check_dtype=False)


@pytest.mark.real_data
def test_known_qc(split_sample_concordance_tables, real_data_cache):
    exp_df = pd.read_csv(
        real_data_cache / "production_outputs/concordance/InternalQcKnown.csv"
    ).set_index(["Sample_ID1", "Sample_ID2"])
    obs_df = (
        split_sample_concordance.read_known_sample_concordance(
            split_sample_concordance_tables / "InternalQcKnown.csv"
        )
        .set_index(["Sample_ID1", "Sample_ID2"])
        .rename(
            {
                "PLINK_PI_HAT": "PI_HAT",
                "PLINK_concordance": "Concordance",
                "Subject_ID1": "Subject_ID",
            },
            axis=1,
        )
        .reindex(exp_df.columns, axis=1)
    )

    assert_frame_equal(exp_df, obs_df, check_like=True, check_dtype=False)


@pytest.mark.real_data
def test_known_study(split_sample_concordance_tables, real_data_cache):
    exp_df = pd.read_csv(
        real_data_cache / "production_outputs/concordance/StudySampleKnown.csv"
    ).set_index(["Sample_ID1", "Sample_ID2"])
    obs_df = (
        split_sample_concordance.read_known_sample_concordance(
            split_sample_concordance_tables / "StudySampleKnown.csv"
        )
        .set_index(["Sample_ID1", "Sample_ID2"])
        .rename(
            {
                "PLINK_PI_HAT": "PI_HAT",
                "PLINK_concordance": "Concordance",
                "Subject_ID1": "Subject_ID",
            },
            axis=1,
        )
        .reindex(exp_df.columns, axis=1)
    )

    assert_frame_equal(exp_df, obs_df, check_like=True, check_dtype=False)


@pytest.mark.real_data
def test_unknown(split_sample_concordance_tables, real_data_cache):
    exp_df = pd.read_csv(
        real_data_cache / "production_outputs/concordance/UnknownReplicates.csv"
    ).set_index(["Sample_ID1", "Sample_ID2"])
    obs_df = (
        split_sample_concordance.read_known_sample_concordance(
            split_sample_concordance_tables / "UnknownReplicates.csv"
        )
        .set_index(["Sample_ID1", "Sample_ID2"])
        .rename({"PLINK_PI_HAT": "PI_HAT", "PLINK_concordance": "Concordance"}, axis=1)
        .reindex(exp_df.columns, axis=1)
    )

    assert_frame_equal(exp_df, obs_df, check_like=True, check_dtype=False)
