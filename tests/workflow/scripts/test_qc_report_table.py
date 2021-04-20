import shutil
from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from cgr_gwas_qc.testing.data import RealData
from cgr_gwas_qc.workflow.scripts import qc_report_table


@pytest.mark.real_data
@pytest.fixture(scope="module")
def qc_report_xlsx(tmp_path_factory) -> Path:
    tmp_path = tmp_path_factory.mktemp("excel_table")
    shutil.copy(
        RealData()
        / "production_outputs/word_doc/SR0446-001_12_QC_Report_1011201995419_casecontrol_20191011.xls",
        tmp_path / "report.xlsx",
    )
    return tmp_path / "report.xlsx"


@pytest.mark.real_data
def test_all_qc(real_sample_sheet_csv, sample_qc_csv, qc_report_xlsx):
    dtypes = {
        "Contaminated": bool,
        "Low Call Rate": bool,
        "Sex Discordant": bool,
        "Unexpected Replicate": bool,
    }

    exp_df = (
        pd.read_excel(qc_report_xlsx, sheet_name="ALL_QC", engine="openpyxl")
        .set_index("Sample_ID")
        .sort_index()
        .reindex(dtypes.keys(), axis=1)
        .fillna(False)
        .astype(dtypes)
    )

    obs_df = (
        qc_report_table._all_qc(real_sample_sheet_csv, sample_qc_csv)
        .set_index("Sample_ID")
        .sort_index()
        .reindex(dtypes.keys(), axis=1)
        .fillna(False)
        .astype(dtypes)
    )

    assert_frame_equal(exp_df, obs_df)


@pytest.fixture
def graf_text(tmp_path):
    outfile = tmp_path / "graf.txt"
    outfile.write_text(
        "DS No.\tSample\t#SNPs\tGD1\tGD2\tGD3\tGD4\tF(%)\tE(%)\tA(%)\tAfrican\tEuropean\tAsian\tMexican\tIndian-Pakistani\n"
        "1\tSB0001\t3366\t1.611165\t1.446184\t0.831981\t0.024308\t3.65\t96.35\t0.00\t1.110809\t0.895608\t1.182723\t0.962664\t0.938356\n"
        "1\tSB0002\t3369\t1.359497\t1.204182\t0.792694\t0.004533\t84.09\t13.00\t2.91\t0.853367\t1.113336\t1.247473\t1.100141\t1.095608\n"
        "1\tSB0003\t3371\t1.618876\t1.447833\t0.836734\t0.037455\t2.14\t97.86\t0.00\t1.119291\t0.895800\t1.179998\t0.971486\t0.934031\n"
        "1\tSB0004\t3371\t1.653709\t1.444371\t0.830352\t0.020863\t0.00\t94.62\t5.38\t1.134156\t0.891016\t1.149434\t0.949329\t0.928466\n"
    )
    return outfile


@pytest.mark.real_data
def test_ancestry(sample_qc_df, graf_text, tmp_path):
    fake_qc = sample_qc_df.head(4)
    fake_qc.Sample_ID = ["SB0001", "SB0002", "SB0003", "SB0004"]
    fake_qc.to_csv(tmp_path / "fake.csv", index=False)
    obs_df = qc_report_table._ancestry(tmp_path / "fake.csv", graf_text)

    assert (4, 16) == obs_df.shape


@pytest.mark.real_data
def test_sample_concordance(sample_qc_csv, sample_concordance_csv):
    obs_df = qc_report_table._sample_concordance(sample_qc_csv, sample_concordance_csv)
    assert qc_report_table._SAMPLE_CONCORDANCE_COLUMNS == obs_df.columns.tolist()
    assert obs_df.iloc[0, :].notna().all()


@pytest.mark.real_data
def test_population_concordance(agg_population_concordance_csv, tmp_path):
    test_file = tmp_path / "test.xlsx"
    with pd.ExcelWriter(test_file) as writer:
        qc_report_table._population_concordance(agg_population_concordance_csv, writer)

    obs_df = pd.read_excel(test_file, "EUR_IBD", engine="openpyxl")
    assert qc_report_table._POPULATION_CONCORDANCE_COLUMNS == obs_df.columns.tolist()
    assert obs_df.iloc[0, :].notna().all()
    assert obs_df.shape[0] == 273


@pytest.fixture
def pop_qc_w_relatives(population_qc_df, tmp_path):
    fake_pop = population_qc_df.copy()

    force_relatives = fake_pop.sample(n=4).Subject_ID
    relative_ids = force_relatives.sort_values()
    mask = fake_pop.Subject_ID.isin(relative_ids)

    fake_pop.loc[mask, "relatives"] = "|".join(relative_ids)
    fake_pop.loc[mask, "QC_Family_ID"] = "fam01"

    outfile = tmp_path / "fake_population_qc.csv"
    fake_pop.to_csv(outfile)
    return outfile
