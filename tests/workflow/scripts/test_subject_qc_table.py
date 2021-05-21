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
