from io import StringIO
from textwrap import dedent

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_series_equal

from cgr_gwas_qc.testing.data import RealData


@pytest.mark.real_data
@pytest.fixture(scope="module")
def ss_df(real_cfg):
    return real_cfg.ss.set_index("Sample_ID")


@pytest.mark.real_data
def test_read_imiss_start(ss_df):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_imiss

    # GIVEN: the call rates from the initial starting point and a list of Sample IDs
    filename = RealData() / "production_outputs/plink_start/samples_start.imiss"

    # WHEN: I parse the imiss table.
    sr = _read_imiss(filename, ss_df.index, "Call_Rate_Initial")

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "Call_Rate_Initial"


@pytest.mark.real_data
def test_read_imiss_cr1(ss_df):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_imiss

    # GIVEN: the call rates after CR1 filters and a list of Sample IDs
    filename = RealData() / "production_outputs/plink_filter_call_rate_1/samples_filter1.imiss"

    # WHEN: I parse the imiss table.
    sr = _read_imiss(filename, ss_df.index, "Call_Rate_1")

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "Call_Rate_1"


@pytest.mark.real_data
def test_read_imiss_cr2(ss_df):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_imiss

    # GIVEN: the call rates after CR2 filters and a list of Sample IDs
    filename = RealData() / "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss"

    # WHEN: I parse the imiss table.
    sr = _read_imiss(filename, ss_df.index, "Call_Rate_2")

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "Call_Rate_2"


@pytest.mark.real_data
def test_read_sexcheck_cr1(ss_df):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_sexcheck_cr1

    # GIVEN: A test sample sheet, plink sexcheck file, and the expected sex
    # calls from the sample table
    filename = RealData() / "production_outputs/plink_filter_call_rate_1/samples_filter1.sexcheck"
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
def test_read_concordance(ss_df, sample_concordance_csv):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_concordance

    # GIVEN: A test sample sheet, config, real production outputs, and a list of Sample_IDs
    # WHEN: Check for discordant replicates
    df = _read_concordance(sample_concordance_csv, ss_df.index)

    # THEN: Basic properties
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "Sample_ID"


@pytest.mark.real_data
@pytest.fixture
def updated_contam(software_params, tmp_path):
    filename = RealData() / "production_outputs/all_contam/contam.csv"
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
def test_read_intensity(ss_df):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_intensity

    # GIVEN: A test sample sheet, real production outputs, and a list of Sample_IDs
    filename = RealData() / "production_outputs/all_sample_idat_intensity/idat_intensity.csv"

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
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _identifiler_reason

    df = pd.DataFrame(
        {
            "identifiler_needed": [False, True, True, True],
            "one": [False, True, False, True],
            "two": [False, True, True, True],
            "three": [False, False, False, True],
        }
    )

    obs_ = _identifiler_reason(df, ["one", "two", "three"])
    exp_ = pd.Series(["", "one;two", "two", "one;two;three"])
    assert_series_equal(obs_, exp_)


def test_find_study_subject_representative():
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _find_study_subject_representative

    # GIVEN: A slimmed down QC report
    df = pd.read_csv(
        StringIO(
            dedent(
                """
        Group_By_Subject_ID,Sample_ID,Call_Rate_2,is_pass_sample_qc
        SP00001,S001,.98,True
        SP00002,S002,.98,True
        SP00002,S003,.99,True
        IB00001,I001,.99,False
        SP00003,S004,.90,False
        SP00003,S005,.91,True
        SP00003,S006,.89,True
        SP00004,S007,.89,False
        SP00004,S008,.9,False
        """
            )
        ),
        index_col="Sample_ID",
    )

    # WHEN: I look for the representative sample for each Subject Group
    # THEN: observed should be the same as expected.
    obs_ = _find_study_subject_representative(df)
    exp_ = pd.Series(
        [True, False, True, False, False, True, False, False, False],
        index=["S001", "S002", "S003", "I001", "S004", "S005", "S006", "S007", "S008"],
    )

    # THEN: Basic properties
    assert_series_equal(obs_, exp_, check_names=False)


def test_find_study_subject_with_no_representative():
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import (
        _find_study_subject_with_no_representative,
    )

    # GIVEN: A slimmed down QC report
    df = pd.read_csv(
        StringIO(
            dedent(
                """
        Group_By_Subject_ID,Sample_ID,is_subject_representative,is_internal_control
        SP00001,S001,True,False
        SP00001,S002,False,False
        SP00002,S003,False,False
        SP00002,S004,False,False
        SP00002,S005,True,False
        SP00003,S006,False,False
        SP00003,S007,False,False
        SP00004,I001,False,True
        """
            )
        ),
        index_col="Sample_ID",
    )

    # WHEN: I look for the representative sample for each Subject Group
    # THEN: observed should be the same as expected.
    obs_ = _find_study_subject_with_no_representative(df)
    exp_ = pd.Series(
        data=[False, False, False, False, False, True, True, False],
        index=["S001", "S002", "S003", "S004", "S005", "S006", "S007", "I001"],
    )

    # THEN: Basic properties
    assert_series_equal(obs_, exp_, check_names=False)
