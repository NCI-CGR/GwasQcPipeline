import pandas as pd
import pytest
from pytest_mock import MockerFixture

from cgr_gwas_qc.reporting.qc_exclusions import ExclusionTables


def test_pre_qc(fake_cfg, mocker: MockerFixture):
    config = mocker.MagicMock()
    config.Sample_IDs_to_remove = ["SB00001_PB0001_A01"]
    pre_qc = ExclusionTables._pre_qc(config, fake_cfg.ss)
    assert 1 == pre_qc.is_array_processing_failure.sum()


def test_pop_qc(fake_cfg):
    df = fake_cfg.ss.assign(Subject_ID=lambda x: x.LIMS_Individual_ID).assign(
        is_extreme_autosomal_heterozygosity=False
    )
    pop_qc = ExclusionTables._pop_qc(df)
    assert (6, 2) == pop_qc.dropna().shape


@pytest.mark.real_data
@pytest.fixture
def exclusion_df(real_cfg, sample_qc_df, population_qc_df) -> pd.DataFrame:
    pre_qc = ExclusionTables._pre_qc(real_cfg.config, real_cfg.ss)
    pop_qc = ExclusionTables._pop_qc(population_qc_df)
    return sample_qc_df.merge(pre_qc, how="outer").merge(pop_qc, how="left")


@pytest.mark.real_data
def test_sample_exclusion_table(exclusion_df):
    expected = (
        "| Filter Reason/Description      |   Controls |   Cases |   QC |   Other |   All Samples |\n"
        "|:-------------------------------|-----------:|--------:|-----:|--------:|--------------:|\n"
        "| Array Processing Failure       |          0 |       0 |    0 |       0 |             0 |\n"
        "| Completion Rate 1st Filter     |          2 |       2 |    0 |       0 |             4 |\n"
        "| Completion Rate 2nd Filter     |          7 |       3 |    0 |       0 |            10 |\n"
        "| Contaminated                   |          0 |       0 |    0 |       0 |             0 |\n"
        "| Internal QC Samples Removed    |          0 |       0 |    7 |       0 |             7 |\n"
        "| Samples Remaining for Analysis |        137 |      62 |    0 |       0 |           199 |\n"
        "| Expected Duplicates Removed    |          8 |       7 |    0 |       0 |            15 |\n"
        "| Subjects Remaining             |        129 |      55 |    0 |       0 |           184 |"
    )

    obs = ExclusionTables._sample_exclusion_table(exclusion_df)

    assert expected == obs


@pytest.mark.real_data
def test_subject_exclusion_table(exclusion_df):
    expected = (
        "| Filter Reason         |   Controls |   Cases |   QC |   Other |   All Subjects |\n"
        "|:----------------------|-----------:|--------:|-----:|--------:|---------------:|\n"
        "| Sex Discordant        |          3 |       0 |    0 |       0 |              3 |\n"
        "| Unexpected Replicates |          0 |       0 |    0 |       0 |              0 |\n"
        "| Autosomal Het         |          0 |       0 |    0 |       0 |              0 |"
    )

    obs = ExclusionTables._subject_exclusion_table(exclusion_df)

    assert expected == obs
