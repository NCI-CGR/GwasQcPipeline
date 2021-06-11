import pytest

from cgr_gwas_qc.reporting.qc_exclusions import ExclusionTables


@pytest.mark.real_data
@pytest.fixture
def exclusion_tables(sample_qc_df, subject_qc_df, agg_population_qc_df) -> ExclusionTables:
    return ExclusionTables.construct(sample_qc_df, subject_qc_df, agg_population_qc_df)


@pytest.mark.real_data
def test_sample_exclusion_table(exclusion_tables):
    expected = (
        "| Filter Reason/Description      |   Controls |   Cases |   QC |   Other |   All Samples |\n"
        "|:-------------------------------|-----------:|--------:|-----:|--------:|--------------:|\n"
        "| User Excluded                  |          0 |       0 |    0 |       0 |             0 |\n"
        "| Array Processing Failure       |          0 |       0 |    0 |       0 |             0 |\n"
        "| Completion Rate 1st Filter     |          2 |       2 |    0 |       0 |             4 |\n"
        "| Completion Rate 2nd Filter     |          7 |       3 |    0 |       0 |            10 |\n"
        "| Contaminated                   |          0 |       0 |    0 |       0 |             0 |\n"
        "| Internal QC Samples Removed    |          0 |       0 |    7 |       0 |             7 |\n"
        "| Samples Remaining for Analysis |        137 |      62 |    0 |       0 |           199 |\n"
        "| Expected Duplicates Removed    |          8 |       7 |    0 |       0 |            15 |\n"
        "| Subjects Remaining             |        129 |      55 |    0 |       0 |           184 |"
    )

    assert expected == exclusion_tables.sample_exclusions


@pytest.mark.real_data
def test_subject_exclusion_table(exclusion_tables):
    expected = (
        "| Filter Reason/Description   |   Controls |   Cases |   QC |   Other |   All Subjects |\n"
        "|:----------------------------|-----------:|--------:|-----:|--------:|---------------:|\n"
        "| Sex Discordant              |          3 |       0 |    0 |       0 |              3 |\n"
        "| Unexpected Replicates       |          0 |       0 |    0 |       0 |              0 |\n"
        "| Autosomal Het               |          0 |       0 |    0 |       0 |              0 |"
    )

    assert expected == exclusion_tables.subject_exclusions
