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


def test_add_case_control(fake_cfg, concordance):
    df = concordance.pipe(sample_concordance._add_case_control, fake_cfg.ss)
    assert "case_control1" in df.columns
    assert "case_control2" in df.columns


def test_add_internal_control(fake_cfg, concordance):
    df = concordance.pipe(sample_concordance._add_internal_control, fake_cfg.ss)
    assert "is_internal_control1" in df.columns
    assert "is_internal_control2" in df.columns


def test_add_expected_replicate(fake_cfg, concordance):
    df = concordance.pipe(sample_concordance._add_expected_replicates, fake_cfg.ss)
    assert 1 == df.is_expected_replicate.sum()


def test_add_expected_replicate_missing_pair(fake_cfg, concordance):
    """Make sure we add missing expected replicates.

    The concordance table does not contain all pairwise combinations for
    storage reasons. I want to make sure if an expected replicate is not in
    the table it gets added.
    """
    # GIVEN: That I have an expected replicate that is not in the concordance table
    ss = fake_cfg.ss.copy()
    ss.loc[0, "replicate_ids"] = "SP00001|SP10000"  # Add an expected replicate
    pseudo_data = pd.DataFrame({"Sample_ID": ["SP10000"], "replicate_ids": ["SP00001|SP10000"]})
    ss = pd.concat([ss, pseudo_data], ignore_index=True)

    # THEN: I should add that pair and flag them as an expected_replicate
    df = concordance.pipe(sample_concordance._add_expected_replicates, ss)
    print("!!!!!!!!")
    print(df["is_expected_replicate"].values)
    print("!!!!!!!!")
    assert 2 == df.is_expected_replicate.sum()


def test_add_discordant_replicate(fake_cfg, concordance):
    df = concordance.pipe(sample_concordance._add_expected_replicates, fake_cfg.ss).pipe(
        sample_concordance._add_discordant_replicates
    )
    assert 0 == df["Expected Replicate Discordance"].sum()


def test_add_discordant_replicate_missing_pair(fake_cfg, concordance):
    """Make sure we add missing expected replicates.

    The concordance table does not contain all pairwise combinations for
    storage reasons. I want to make sure if an expected replicate is not in
    the table it gets added.
    """
    # GIVEN: That I force the replciate to be discordant
    concordance.loc[("SP00002", "SP00003"), "PLINK_is_ge_concordance"] = False

    # THEN: I should have one flag
    df = concordance.pipe(sample_concordance._add_expected_replicates, fake_cfg.ss).pipe(
        sample_concordance._add_discordant_replicates
    )
    assert 1 == df["Expected Replicate Discordance"].sum()


@pytest.mark.parametrize(
    "expected_result,is_rep,plink,graf,king",
    [
        (False, False, True, "ID", "ID"),  # Not an expected replicate
        # PLINK Concordant
        (False, True, True, "ID", "ID"),  # all concordant
        (False, True, True, "UN", "UN"),  # GRAF and KING don't matter
        (False, True, True, pd.NA, pd.NA),  # GRAF and KING don't matter
        # PLINK Discordant
        (True, True, False, "ID", "UN"),  # GRAF and KING don't matter
        (True, True, False, "UN", "ID"),  # GRAF and KING don't matter
        (True, True, False, pd.NA, pd.NA),  # GRAF and KING don't matter
        (True, True, False, "UN", "UN"),  # GRAF and KING don't matter
        # PLINK Missing
        (False, True, pd.NA, "ID", "ID"),  # GRAF and KING concordant
        (False, True, pd.NA, "ID", "UN"),  # GRAF and KING disagree
        (False, True, pd.NA, "UN", "ID"),  # GRAF and KING disagree
        (False, True, pd.NA, "ID", pd.NA),  # GRAF concordant
        (False, True, pd.NA, pd.NA, "ID"),  # KING concordant
        (True, True, pd.NA, "UN", pd.NA),  # GRAF discordant
        (True, True, pd.NA, pd.NA, "UN"),  # KING discordant
        (True, True, pd.NA, pd.NA, pd.NA),  # All missing: flag as discordant replicate
    ],
)
def test_add_discordant_logic(expected_result, is_rep, plink, graf, king):
    df = pd.DataFrame(
        {
            "is_expected_replicate": [is_rep],
            "PLINK_is_ge_concordance": [plink],
            "GRAF_relationship": [graf],
            "KING_relationship": [king],
        }
    )

    result = sample_concordance._add_discordant_replicates(df).is_discordant_replicate.squeeze()
    assert expected_result == result


def test_add_unexpected_replicate_none(fake_cfg, concordance):
    # GIVEN: Test data does not have any unexpected replicates
    df = concordance.pipe(sample_concordance._add_expected_replicates, fake_cfg.ss).pipe(
        sample_concordance._add_unexpected_replicates
    )

    # THEN: I should have all False
    assert 0 == df.is_unexpected_replicate.sum()


def test_add_unexpected_replicate_one(fake_cfg, concordance):
    # GIVEN: I add an unexpected replicate
    concordance.loc[0, "PLINK_is_ge_concordance"] = True
    df = concordance.pipe(sample_concordance._add_expected_replicates, fake_cfg.ss).pipe(
        sample_concordance._add_unexpected_replicates
    )
    # THEN: I should have one sample flagged
    assert 1 == df.is_unexpected_replicate.sum()


@pytest.mark.parametrize(
    "expected_result,is_rep,plink,graf,king",
    [
        (False, True, True, "ID", "ID"),  # expected replicate
        (True, False, True, "ID", "ID"),  # All concordant
        (True, False, True, "UN", "UN"),  # plink concordant
        (True, False, False, "ID", "UN"),  # graf concordant
        (True, False, False, "UN", "ID"),  # king concordant
        (True, False, True, pd.NA, pd.NA),  # plink concordant with missing
        (True, False, pd.NA, "ID", pd.NA),  # graf concordant with missing
        (True, False, pd.NA, pd.NA, "ID"),  # king concordant with missing
        (False, False, False, "UN", "UN"),  # all discordant
        (False, False, pd.NA, pd.NA, pd.NA),  # all missing
    ],
)
def test_add_unexpected_logic(expected_result, is_rep, plink, graf, king):
    df = pd.DataFrame(
        {
            "is_expected_replicate": [is_rep],
            "PLINK_is_ge_concordance": [plink],
            "GRAF_relationship": [graf],
            "KING_relationship": [king],
        }
    )

    result = sample_concordance._add_unexpected_replicates(df).is_unexpected_replicate.squeeze()
    assert expected_result == result
