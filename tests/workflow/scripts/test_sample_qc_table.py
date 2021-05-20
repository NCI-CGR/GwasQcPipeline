from textwrap import dedent

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_series_equal

from cgr_gwas_qc.workflow.scripts import sample_qc_table


@pytest.mark.real_data
@pytest.fixture(scope="module")
def ss_df(real_cfg):
    return real_cfg.ss.set_index("Sample_ID")


@pytest.mark.real_data
def test_read_imiss_start(real_data_cache, ss_df):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_imiss

    # GIVEN: the call rates from the initial starting point and a list of Sample IDs
    filename = real_data_cache / "legacy_outputs/plink_start/samples_start.imiss"

    # WHEN: I parse the imiss table.
    sr = _read_imiss(filename, ss_df.index, "Call_Rate_Initial")

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "Call_Rate_Initial"


@pytest.mark.real_data
def test_read_imiss_cr1(real_data_cache, ss_df):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_imiss

    # GIVEN: the call rates after CR1 filters and a list of Sample IDs
    filename = real_data_cache / "legacy_outputs/plink_filter_call_rate_1/samples_filter1.imiss"

    # WHEN: I parse the imiss table.
    sr = _read_imiss(filename, ss_df.index, "Call_Rate_1")

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "Call_Rate_1"


@pytest.mark.real_data
def test_read_imiss_cr2(real_data_cache, ss_df):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_imiss

    # GIVEN: the call rates after CR2 filters and a list of Sample IDs
    filename = real_data_cache / "legacy_outputs/plink_filter_call_rate_2/samples_filter2.imiss"

    # WHEN: I parse the imiss table.
    sr = _read_imiss(filename, ss_df.index, "Call_Rate_2")

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "Call_Rate_2"


@pytest.mark.real_data
def test_read_sexcheck_cr1(real_data_cache, ss_df):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_sexcheck_cr1

    # GIVEN: A test sample sheet, plink sexcheck file, and the expected sex
    # calls from the sample table
    filename = real_data_cache / "legacy_outputs/plink_filter_call_rate_1/samples_filter1.sexcheck"
    expected_sex_calls = ss_df["expected_sex"]

    # add an extra sample to check how missing values are handled
    expected_sex_calls["SB_missing_sample"] = "F"

    # WHEN: read the sexcheck table
    df = _read_sexcheck_cr1(filename, expected_sex_calls)

    # THEN: Basic properties
    assert df.index.name == "Sample_ID"
    assert df.predicted_sex.dtype == pd.CategoricalDtype(categories=["M", "F", "U"])
    assert isinstance(df.is_sex_discordant.dtype, pd.BooleanDtype)
    assert df.X_inbreeding_coefficient.dtype is np.dtype("float64")


@pytest.mark.real_data
def test_read_ancestry_GRAF(tmp_path):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_ancestry

    # GIVEN: A test sample sheet, example outputs from GRAF -pop, and a list of Sample_IDs
    filename = tmp_path / "graf.tsv"
    filename.write_text(
        dedent(
            """
            Subject\tSelf-reported ancestry\tGD1\tGD2\tGD3\tGD4\tP_f (%)\tP_e (%)\tP_a (%)\tPopID\tComputed population
            SB034307_PC26469_B09\t\t1.63\t1.45\t0.81\t0.03\t0.6\t99.1\t0.2\t1\tEuropean
            SB034327_PC26469_C09\t\t1.64\t1.47\t0.85\t0.02\t0.0\t100.0\t0.0\t1\tEuropean
            """
        )
    )
    sample_ids = pd.Index(["SB034307_PC26469_B09", "SB034327_PC26469_C09"], name="Sample_ID")

    # WHEN: Parse the GRAF output
    df = _read_ancestry(filename, sample_ids)

    # THEN: Basic properties
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "Sample_ID"
    assert "AFR" in df.columns
    assert "EUR" in df.columns
    assert "ASN" in df.columns
    assert "Ancestry" in df.columns
    assert df.dtypes["AFR"] is np.dtype("float64")
    assert df.dtypes["EUR"] is np.dtype("float64")
    assert df.dtypes["ASN"] is np.dtype("float64")

    # Check the conversion from a percent to a proportion
    assert df.loc["SB034307_PC26469_B09", "EUR"] == 0.991
    assert df.loc["SB034327_PC26469_C09", "EUR"] == 1.0


@pytest.mark.real_data
def test_read_concordance(real_data_cache, ss_df):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_concordance

    # GIVEN: A test sample sheet, config, real production outputs, and a list of Sample_IDs
    # WHEN: Check for discordant replicates
    sample_concordance_csv = real_data_cache / "dev_outputs/sample_level/concordance/summary.csv"
    df = _read_concordance(sample_concordance_csv, ss_df.index)

    # THEN: Basic properties
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "Sample_ID"


@pytest.mark.real_data
@pytest.fixture
def updated_contam(real_data_cache, software_params, tmp_path):
    filename = real_data_cache / "legacy_outputs/all_contam/contam.csv"
    outfile = tmp_path / "contam.csv"
    (
        pd.read_csv(filename)
        .rename({"ID": "Sample_ID"}, axis=1)
        .assign(is_contaminated=lambda x: x["%Mix"] > software_params.contam_threshold)
        .to_csv(outfile, index=False)
    )
    return outfile


@pytest.mark.real_data
def test_read_contam(ss_df, updated_contam):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_contam

    # GIVEN: A test sample sheet, config, real production outputs, and a list of Sample_IDs
    # WHEN: Parse the contamination table
    df = _read_contam(updated_contam, ss_df.index)

    # THEN: Basic properties
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "Sample_ID"
    assert "Contamination_Rate" in df.columns
    assert "is_contaminated" in df.columns


@pytest.mark.real_data
def test_read_contam_file_name_none(ss_df):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_contam

    # GIVEN: A test sample sheet, config, real production outputs, and a list of Sample_IDs
    filename = None

    # WHEN: Parse the contamination table.
    df = _read_contam(filename, ss_df.index)

    # THEN: Basic properties
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "Sample_ID"
    assert "Contamination_Rate" in df.columns
    assert "is_contaminated" in df.columns


@pytest.mark.real_data
def test_read_intensity(real_data_cache, ss_df):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_intensity

    # GIVEN: A test sample sheet, real production outputs, and a list of Sample_IDs
    filename = real_data_cache / "legacy_outputs/all_sample_idat_intensity/idat_intensity.csv"

    # WHEN: Parse the read intensity table.
    sr = _read_intensity(filename, ss_df.index)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "IdatIntensity"
    assert sr.dtype is np.dtype("float64")


@pytest.mark.real_data
def test_read_intensity_file_name_none(ss_df):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_intensity

    # GIVEN: A test sample sheet, No production outputs, and a list of Sample_IDs
    filename = None

    # WHEN: Parse the read intensity table with file name of None
    sr = _read_intensity(filename, ss_df.index)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "IdatIntensity"


def test_identifiler_reason():
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _get_reason

    df = pd.DataFrame(
        {
            "identifiler_needed": [False, True, True, True],
            "one": [False, True, False, True],
            "two": [False, True, True, True],
            "three": [False, False, False, True],
        }
    )
    flags = {"one": "A", "two": "B", "three": "C"}
    obs_ = _get_reason(df, flags)
    exp_ = pd.Series([pd.NA, "A|B", "B", "A|B|C"])
    assert_series_equal(obs_, exp_)


@pytest.fixture
def fake_sample_qc() -> pd.DataFrame:
    columns = [
        "Sample_ID",
        "Group_By_Subject_ID",
        "is_sample_exclusion",
        "is_internal_control",
        "Call_Rate_2",
        "is_cr1_filtered",
        "is_cr2_filtered",
        "is_contaminated",
        "is_discordant_replicate",
    ]
    data = [
        ("SP00001", "SB00001", False, False, 0.99, False, False, False, False),
        ("SP00002", "SB00002", False, False, 0.82, False, True, False, False),
        ("SP00003", "SB00003", False, False, 0.99, False, False, True, True),
        ("SP00004", "SB00003", False, False, 0.99, False, False, False, True),
        ("SP00005", "SB00004", False, False, 0.99, False, False, False, True),
        ("SP00006", "SB00004", False, False, 0.99, False, False, False, True),
        ("SP00007", "SB00005", False, False, 0.99, False, False, False, False),
        ("SP00008", "SB00006", False, False, 0.99, False, False, False, False),
        ("SP00009", "SB00007", False, False, 0.99, False, False, False, False),
        ("SP00010", "SB00008", False, False, 0.99, False, False, False, False),
        ("SP00011", "SB00008", False, False, 0.94, False, False, False, False),
        ("SP00012", "SB00009", False, False, 0.99, False, False, False, False),
        ("SP00013", "SB00009", False, False, 0.98, False, False, False, False),
        ("SP00014", "SB00009", False, False, 0.97, False, False, False, False),
    ]
    return pd.DataFrame(data, columns=columns).set_index("Sample_ID")


@pytest.mark.parametrize(
    "contam,rep_discordant,num_removed",
    [(False, False, 1), (True, False, 2), (False, True, 5), (True, True, 5)],  # call rate filtered
)
def test_add_analytic_exclusion(fake_sample_qc, contam, rep_discordant, num_removed):
    sample_qc_table._add_analytic_exclusion(fake_sample_qc, contam, rep_discordant)
    assert num_removed == fake_sample_qc.analytic_exclusion.sum()


@pytest.mark.parametrize(
    "contam,rep_discordant,num_subjects",
    [(False, False, 8), (True, False, 8), (False, True, 6), (True, True, 6)],
)
def test_add_subject_representative(fake_sample_qc, contam, rep_discordant, num_subjects):
    sample_qc_table._add_analytic_exclusion(fake_sample_qc, contam, rep_discordant)
    sample_qc_table._add_subject_representative(fake_sample_qc)
    assert num_subjects == fake_sample_qc.is_subject_representative.sum()


@pytest.mark.parametrize(
    "contam,rep_discordant,num_subjects",
    [(False, False, 1), (True, False, 1), (False, True, 3), (True, True, 3)],
)
def test_add_subject_dropped_from_study(fake_sample_qc, contam, rep_discordant, num_subjects):
    sample_qc_table._add_analytic_exclusion(fake_sample_qc, contam, rep_discordant)
    sample_qc_table._add_subject_representative(fake_sample_qc)
    sample_qc_table._add_subject_dropped_from_study(fake_sample_qc)
    subs = fake_sample_qc.groupby("Group_By_Subject_ID").subject_dropped_from_study.all()
    assert num_subjects == subs.sum()
