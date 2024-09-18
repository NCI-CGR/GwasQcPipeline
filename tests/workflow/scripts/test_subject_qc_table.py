import pandas as pd
import pytest

from cgr_gwas_qc.workflow.scripts import subject_qc_table


@pytest.mark.real_data
def test_subject_qc_table(real_data_cache, sample_qc_df):
    sample_concordance_csv = real_data_cache / "dev_outputs/sample_level/concordance/summary.csv"
    df = subject_qc_table._sample_qc_to_subject_qc(sample_qc_df).pipe(
        subject_qc_table._add_unexpected_replicate_ids, sample_concordance_csv
    )
    assert 0 == df.is_unexpected_replicate.sum()
    assert 3 == df.is_sex_discordant.sum()


@pytest.fixture
def fake_sample_qc() -> pd.DataFrame:
    columns = [
        "Sample_ID",
        "Group_By_Subject_ID",
        "is_contaminated",
        "is_unexpected_replicate",
        "unexpected_replicate_ids",
    ]
    data = [
        ("SP00001", "SB00001", True, True, "SB00001|SB00002"),
        ("SP00002", "SB00002", False, True, "SB00001|SB00002"),
        ("SP00003", "SB00003", False, False, ""),
        ("SP00004", "SB00004", True, False, ""),
        ("SP00005", "SB00005", False, True, "SB00005|SB00006"),
        ("SP00006", "SB00006", False, True, "SB00005|SB00006"),
    ]
    return pd.DataFrame(data, columns=columns)


# create a test case that tests the new column "unexpected_replicate_status"
@pytest.mark.parametrize(
    "subject_id, predicted_status, predicted_unexp_rep",
    [
        ("SB00001", 2, True),
        ("SB00002", 1, False),
        ("SB00003", 0, False),
        ("SB00004", 0, False),
        ("SB00005", 3, True),
        ("SB00006", 3, True),
    ],
)
def test_add_unexpected_replicate_status(
    fake_sample_qc, subject_id, predicted_status, predicted_unexp_rep
):
    df = subject_qc_table._add_unexpected_replicate_status(fake_sample_qc).copy()
    df = df.set_index("Group_By_Subject_ID")

    assert (
        df.loc[subject_id, "unexpected_replicate_status"] == predicted_status
        and df.loc[subject_id, "is_unexpected_replicate"] == predicted_unexp_rep
    )
