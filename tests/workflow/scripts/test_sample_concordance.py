import pandas as pd
import pytest

from cgr_gwas_qc.testing.data import FakeData
from cgr_gwas_qc.workflow.scripts import sample_concordance


@pytest.fixture
def concordance():
    data_cache = FakeData()
    return sample_concordance.build(
        data_cache / "cgr/concordance.csv",
        data_cache / "graf/relatedness.tsv",
        data_cache / "king/king.kin0",
    )


def test_build(concordance):
    assert 4 == concordance.filter(regex="^PLINK").shape[1]
    assert 3 == concordance.filter(regex="^GRAF").shape[1]
    assert 2 == concordance.filter(regex="^KING").shape[1]
    assert ["Sample_ID1", "Sample_ID2"] == concordance.index.names


def test_add_subject(fake_cfg, concordance):
    df = concordance.pipe(sample_concordance._add_subject, fake_cfg.ss)
    assert "Subject_ID1" in df.columns
    assert "Subject_ID2" in df.columns


def test_add_internal_control(fake_cfg, concordance):
    df = concordance.pipe(sample_concordance._add_internal_control, fake_cfg.ss)
    assert "is_internal_control1" in df.columns
    assert "is_internal_control2" in df.columns


def test_add_expected_replicate(fake_cfg, concordance):
    df = concordance.pipe(sample_concordance._add_expected_replicates, fake_cfg.ss)
    assert 1 == df.expected_replicate.sum()


def test_add_expected_replicate_missing_pair(fake_cfg, concordance):
    """Make sure we add missing expected replicates.

    The concordance table does not contain all pairwise combinations for
    storage reasons. I want to make sure if an expected replicate is not in
    the table it gets added.
    """
    # GIVEN: That I have an expected replicate that is not in the concordance table
    ss = fake_cfg.ss.copy()
    ss.loc[0, "replicate_ids"] = "SP00001|SP10000"  # Add an expected replicate
    ss = ss.append(
        pd.Series({"Sample_ID": "SP10000", "replicate_ids": "SP00001|SP10000"}), ignore_index=True
    )

    # THEN: I should add that pair and flag them as an expected_replicate
    df = concordance.pipe(sample_concordance._add_expected_replicates, ss)
    assert 2 == df.expected_replicate.sum()


def test_add_unexpected_replicate_none(fake_cfg, concordance):
    # GIVEN: Test data does not have any unexpected replicates
    df = concordance.pipe(sample_concordance._add_expected_replicates, fake_cfg.ss).pipe(
        sample_concordance._add_unexpected_replicates
    )

    # THEN: I should have all False
    assert 0 == df.unexpected_replicate.sum()


def test_add_unexpected_replicate_one(fake_cfg, concordance):
    # GIVEN: I add an unexpected replicate
    concordance.loc[0, "PLINK_is_ge_concordance"] = True
    df = concordance.pipe(sample_concordance._add_expected_replicates, fake_cfg.ss).pipe(
        sample_concordance._add_unexpected_replicates
    )
    # THEN: I should have one sample flagged
    assert 1 == df.unexpected_replicate.sum()
