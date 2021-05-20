import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from cgr_gwas_qc.workflow.scripts import split_sample_concordance


@pytest.mark.real_data
def test_known_full(real_data_cache):
    legacy_file = real_data_cache / "legacy_outputs/concordance/KnownReplicates.csv"
    dev_file = real_data_cache / "dev_outputs/sample_level/concordance/KnownReplicates.csv"
    legacy_df = pd.read_csv(legacy_file).set_index(["Sample_ID1", "Sample_ID2"])

    dev_df = (
        split_sample_concordance.read_known_sample_concordance(dev_file)
        .set_index(["Sample_ID1", "Sample_ID2"])
        .rename(
            {
                "PLINK_PI_HAT": "PI_HAT",
                "PLINK_concordance": "Concordance",
                "Subject_ID1": "Subject_ID",
            },
            axis=1,
        )
        .reindex(legacy_df.columns, axis=1)
    )

    assert_frame_equal(legacy_df, dev_df, check_like=True, check_dtype=False)


@pytest.mark.real_data
def test_known_qc(real_data_cache):
    legacy_file = real_data_cache / "legacy_outputs/concordance/InternalQcKnown.csv"
    dev_file = real_data_cache / "dev_outputs/sample_level/concordance/InternalQcKnown.csv"

    legacy_df = pd.read_csv(legacy_file).set_index(["Sample_ID1", "Sample_ID2"])
    dev_df = (
        split_sample_concordance.read_known_sample_concordance(dev_file)
        .set_index(["Sample_ID1", "Sample_ID2"])
        .rename(
            {
                "PLINK_PI_HAT": "PI_HAT",
                "PLINK_concordance": "Concordance",
                "Subject_ID1": "Subject_ID",
            },
            axis=1,
        )
        .reindex(legacy_df.columns, axis=1)
    )

    assert_frame_equal(legacy_df, dev_df, check_like=True, check_dtype=False)


@pytest.mark.real_data
def test_known_study(real_data_cache):
    legacy_file = real_data_cache / "legacy_outputs/concordance/StudySampleKnown.csv"
    dev_file = real_data_cache / "dev_outputs/sample_level/concordance/StudySampleKnown.csv"

    legacy_df = pd.read_csv(legacy_file).set_index(["Sample_ID1", "Sample_ID2"])
    dev_df = (
        split_sample_concordance.read_known_sample_concordance(dev_file)
        .set_index(["Sample_ID1", "Sample_ID2"])
        .rename(
            {
                "PLINK_PI_HAT": "PI_HAT",
                "PLINK_concordance": "Concordance",
                "Subject_ID1": "Subject_ID",
            },
            axis=1,
        )
        .reindex(legacy_df.columns, axis=1)
    )

    assert_frame_equal(legacy_df, dev_df, check_like=True, check_dtype=False)


@pytest.mark.real_data
def test_unknown(real_data_cache):
    legacy_file = real_data_cache / "legacy_outputs/concordance/UnknownReplicates.csv"
    dev_file = real_data_cache / "dev_outputs/sample_level/concordance/UnknownReplicates.csv"

    legacy_df = pd.read_csv(legacy_file).set_index(["Sample_ID1", "Sample_ID2"])
    dev_df = (
        split_sample_concordance.read_known_sample_concordance(dev_file)
        .set_index(["Sample_ID1", "Sample_ID2"])
        .rename({"PLINK_PI_HAT": "PI_HAT", "PLINK_concordance": "Concordance"}, axis=1)
        .reindex(legacy_df.columns, axis=1)
    )

    assert_frame_equal(legacy_df, dev_df, check_like=True, check_dtype=False)
