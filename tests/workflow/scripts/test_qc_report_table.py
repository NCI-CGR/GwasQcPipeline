import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from cgr_gwas_qc.workflow.scripts import qc_report_table


@pytest.mark.real_data
def test_sample_qc(real_data_cache):
    dtypes = {
        "Contaminated": bool,
        "Low Call Rate": bool,
        "Sex Discordant": bool,
        "Unexpected Replicate": bool,
    }

    legacy_file = (
        real_data_cache
        / "legacy_outputs/word_doc/SR0446-001_12_QC_Report_1011201995419_casecontrol_20191011.xlsx"
    )

    legacy_df = (
        pd.read_excel(legacy_file, sheet_name="ALL_QC", engine="openpyxl")
        .set_index("Sample_ID")
        .sort_index()
        .reindex(dtypes.keys(), axis=1)
        .fillna(False)
        .astype(dtypes)
    )

    dev_df = (
        qc_report_table._sample_qc(
            real_data_cache / "dev_outputs/cgr_sample_sheet.csv",
            real_data_cache / "dev_outputs/sample_level/sample_qc.csv",
        )
        .set_index("Sample_ID")
        .sort_index()
        .reindex(dtypes.keys(), axis=1)
        .fillna(False)
        .astype(dtypes)
    )

    assert_frame_equal(legacy_df, dev_df)


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

    assert (4, 18) == obs_df.shape


@pytest.mark.real_data
def test_sample_concordance(real_data_cache):
    sample_qc_csv = real_data_cache / "dev_outputs/sample_level/sample_qc.csv"
    sample_concordance_csv = real_data_cache / "dev_outputs/sample_level/concordance/summary.csv"
    obs_df = qc_report_table._sample_concordance(sample_qc_csv, sample_concordance_csv)
    assert qc_report_table._SAMPLE_CONCORDANCE_COLUMNS == obs_df.columns.tolist()


@pytest.mark.real_data
def test_population_concordance(real_data_cache, tmp_path):
    agg_population_concordance_csv = real_data_cache / "dev_outputs/subject_level/concordance.csv"
    test_file = tmp_path / "test.xlsx"
    with pd.ExcelWriter(test_file) as writer:
        qc_report_table._population_concordance(agg_population_concordance_csv, writer)

    obs_df = pd.read_excel(test_file, "European_IBD", engine="openpyxl")
    assert qc_report_table._POPULATION_CONCORDANCE_COLUMNS == obs_df.columns.tolist()
    assert obs_df.shape[0] == 45


@pytest.mark.real_data
def test_pca(real_data_cache, tmp_path):
    filename = real_data_cache / "dev_outputs/subject_level/population_qc.csv"
    test_file = tmp_path / "test.xlsx"
    with pd.ExcelWriter(test_file) as writer:
        qc_report_table._pca(filename, writer)

    obs_df = pd.read_excel(test_file, "European_PCA", engine="openpyxl")
    assert qc_report_table._PCA_COLUMNS == obs_df.columns.tolist()
    assert not obs_df.isna().all(axis=0).any()  # no empty columns


@pytest.mark.real_data
def test_het(real_data_cache, tmp_path):
    filename = real_data_cache / "dev_outputs/subject_level/population_qc.csv"
    test_file = tmp_path / "test.xlsx"
    with pd.ExcelWriter(test_file) as writer:
        qc_report_table._het(filename, writer)

    obs_df = pd.read_excel(test_file, "European_HET", engine="openpyxl")
    assert qc_report_table._HET_COLUMNS == obs_df.columns.tolist()
    assert not obs_df.isna().all(axis=0).any()  # no empty columns
