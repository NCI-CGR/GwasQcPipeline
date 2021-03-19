from io import StringIO
from textwrap import dedent

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_series_equal

from cgr_gwas_qc import load_config
from cgr_gwas_qc.models.config.software_params import SoftwareParams
from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.testing.data import FakeData, RealData
from cgr_gwas_qc.workflow.scripts.sample_qc_table import CASE_CONTROL_DTYPE, SEX_DTYPE


@pytest.fixture
def sample_ids_short(sample_sheet_short) -> pd.Index:
    return pd.Index(sample_sheet_short.Sample_ID)


@pytest.fixture
def sample_ids_full(sample_sheet_full) -> pd.Index:
    return pd.Index(sample_sheet_full.Sample_ID)


@pytest.mark.real_data
@pytest.mark.parametrize("expected_sex_col", ["Expected_Sex", "Identifiler_Sex"])
def test_wrangle_sample_sheet(sample_sheet_full, expected_sex_col):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _wrangle_sample_sheet

    # GIVEN: A sample sheet and a column to use as `expected_sex_col`
    # WHEN: I wrangle the sample sheet
    ss = _wrangle_sample_sheet(sample_sheet_full, expected_sex_col)

    # THEN: Basic properties
    assert ss.index.name == "Sample_ID"
    assert isinstance(ss.is_internal_control.dtype, pd.BooleanDtype)
    assert ss.expected_sex.dtype == SEX_DTYPE
    assert ss.case_control.dtype == CASE_CONTROL_DTYPE

    # The `expected_sex` column should be the same as the column passed as `expected_sex_col`
    control_Sample_ID = ss.query("case_control == 'QC'").index
    assert_series_equal(
        ss.expected_sex.drop(control_Sample_ID).sort_index().astype("object"),
        sample_sheet_full.set_index("Sample_ID")[expected_sex_col]
        .drop(control_Sample_ID)
        .sort_index(),
        check_names=False,
    )


@pytest.mark.real_data
def test_check_preflight_none(sample_ids_full):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _check_preflight

    sr = _check_preflight(None, sample_ids_full)
    assert not any(sr)


@pytest.mark.real_data
def test_check_preflight(sample_ids_full):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _check_preflight

    sr = _check_preflight(sample_ids_full[:2], sample_ids_full)

    assert "Sample_ID" == sr.index.name
    assert "is_preflight_exclusion" == sr.name
    assert 2 == sum(sr)


@pytest.mark.real_data
def test_read_imiss_start(sample_ids_full):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_imiss

    # GIVEN: the call rates from the initial starting point and a list of Sample IDs
    filename = RealData() / "production_outputs/plink_start/samples_start.imiss"

    # WHEN: I parse the imiss table.
    sr = _read_imiss(filename, sample_ids_full, "Call_Rate_Initial")

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "Call_Rate_Initial"


@pytest.mark.real_data
def test_read_imiss_cr1(sample_ids_full):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_imiss

    # GIVEN: the call rates after CR1 filters and a list of Sample IDs
    filename = RealData() / "production_outputs/plink_filter_call_rate_1/samples_filter1.imiss"

    # WHEN: I parse the imiss table.
    sr = _read_imiss(filename, sample_ids_full, "Call_Rate_1")

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "Call_Rate_1"


@pytest.mark.real_data
def test_read_imiss_cr2(sample_ids_full):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_imiss

    # GIVEN: the call rates after CR2 filters and a list of Sample IDs
    filename = RealData() / "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss"

    # WHEN: I parse the imiss table.
    sr = _read_imiss(filename, sample_ids_full, "Call_Rate_2")

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "Call_Rate_2"


@pytest.mark.real_data
def test_read_sexcheck_cr1(sample_sheet_full):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_sexcheck_cr1

    # GIVEN: A test sample sheet, plink sexcheck file, and the expected sex
    # calls from the sample table
    filename = RealData() / "production_outputs/plink_filter_call_rate_1/samples_filter1.sexcheck"
    expected_sex_calls = sample_sheet_full.set_index("Sample_ID")["Expected_Sex"]
    expected_sex_calls[
        "SB_missing_sample"
    ] = "F"  # add an extra sample to check how missing values are handled

    # WHEN: read the sexcheck table
    df = _read_sexcheck_cr1(filename, expected_sex_calls)

    # THEN: Basic properties
    assert df.index.name == "Sample_ID"
    assert df.predicted_sex.dtype == pd.CategoricalDtype(categories=["M", "F", "U"])
    assert isinstance(df.is_sex_discordant.dtype, pd.BooleanDtype)
    assert df.X_inbreeding_coefficient.dtype is np.dtype("float64")


@pytest.mark.real_data
def test_read_ancestry_GRAF(sample_ids_short, tmp_path):
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

    # WHEN: Parse the GRAF output
    df = _read_ancestry(filename, sample_ids_short)

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
def test_read_known_replicates(sample_ids_full):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_known_replicates

    # GIVEN: A test sample sheet, config, real production outputs, and a list of Sample_IDs
    filename = RealData() / "production_outputs/concordance/KnownReplicates.csv"
    cutoff = SoftwareParams().dup_concordance_cutoff

    # WHEN: Check for discordant replicates
    sr = _read_known_replicates(filename, cutoff, sample_ids_full)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "is_replicate_discordant"
    assert sr.dtype is np.dtype("bool")


@pytest.mark.real_data
def test_read_unknown_replicates(sample_ids_full):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_unknown_replicates

    # GIVEN: A test sample sheet, real production outputs, and a list of Sample_IDs
    filename = RealData() / "production_outputs/concordance/UnknownReplicates.csv"

    # WHEN: I parse unknown concordant samples table.
    sr = _read_unknown_replicates(filename, sample_ids_full)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "is_unexpected_replicate"
    assert sr.dtype is np.dtype("bool")


@pytest.mark.real_data
def test_read_contam(sample_ids_full):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_contam

    # GIVEN: A test sample sheet, config, real production outputs, and a list of Sample_IDs
    filename = RealData() / "production_outputs/all_contam/contam.csv"
    cutoff = SoftwareParams().contam_threshold

    # WHEN: Parse the contamination table
    df = _read_contam(filename, cutoff, sample_ids_full)

    # THEN: Basic properties
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "Sample_ID"
    assert "Contamination_Rate" in df.columns
    assert "is_contaminated" in df.columns
    assert df.dtypes["is_contaminated"] is np.dtype("bool")
    assert df.dtypes["Contamination_Rate"] is np.dtype("float64")


@pytest.mark.real_data
def test_read_contam_file_name_none(sample_ids_full):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_contam

    # GIVEN: A test sample sheet, config, real production outputs, and a list of Sample_IDs
    filename = None
    cutoff = SoftwareParams().contam_threshold

    # WHEN: Parse the contamination table.
    df = _read_contam(filename, cutoff, sample_ids_full)

    # THEN: Basic properties
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "Sample_ID"
    assert "Contamination_Rate" in df.columns
    assert "is_contaminated" in df.columns


@pytest.mark.real_data
def test_read_intensity(sample_ids_full):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_intensity

    # GIVEN: A test sample sheet, real production outputs, and a list of Sample_IDs
    filename = RealData() / "production_outputs/all_sample_idat_intensity/idat_intensity.csv"

    # WHEN: Parse the read intensity table.
    sr = _read_intensity(filename, sample_ids_full)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "IdatIntensity"
    assert sr.dtype is np.dtype("float64")


@pytest.mark.real_data
def test_read_intensity_file_name_none(sample_ids_full):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_intensity

    # GIVEN: A test sample sheet, No production outputs, and a list of Sample_IDs
    filename = None

    # WHEN: Parse the read intensity table with file name of None
    sr = _read_intensity(filename, sample_ids_full)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "IdatIntensity"


def test_check_idat_files(tmp_path):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _check_idats_files

    # GIVEN: Fake data using the gtc entrypoint to create the idat files
    (FakeData(tmp_path).add_user_files("gtc").make_config())

    # WHEN: Scan the folder for Idat files.
    with chdir(tmp_path):
        cfg = load_config()
        sr = _check_idats_files(cfg.ss.set_index("Sample_ID"), cfg.config.user_files.idat_pattern)

    # THEN: Basic properties
    assert sr.index.name == "Sample_ID"
    assert sr.name == "idats_exist"

    # All idat files should be found
    assert all(sr)


@pytest.mark.real_data
def test_check_idat_files_one_missing(tmp_path):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _check_idats_files

    # GIVEN: Fake data using the gtc entrypoint to create the idat files
    (FakeData(tmp_path).add_user_files("gtc").make_config())

    # WHEN: I add an extra sample (without idat file) and scan the folder for Idat files.
    with chdir(tmp_path):
        cfg = load_config()
        fake_record = cfg.ss.iloc[0, :].copy()
        fake_record["Sample_ID"] = "fake_Sample_ID"
        fake_record["SentrixBarcode_A"] = "fake_barcode"
        cfg.ss = cfg.ss.append(fake_record, ignore_index=True)
        sr = _check_idats_files(cfg.ss.set_index("Sample_ID"), cfg.config.user_files.idat_pattern)

    # I should have 4 samples with idat files found and 1 missing these files.
    assert sum(sr) == 4
    assert sum(~sr) == 1


def test_identifiler_reason():
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _identifiler_reason

    df = pd.DataFrame(
        {
            "Identifiler_Needed": [False, True, True, True],
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
        Group_By_Subject_ID,Sample_ID,Call_Rate_2,is_call_rate_filtered,is_contaminated,is_replicate_discordant,is_internal_control
        SB001,S001,.98,False,False,False,False
        SB002,S002,.98,False,False,False,False
        SB002,S003,.99,False,False,False,False
        IB001,I001,.99,False,False,False,True
        SB003,S004,.90,False,True,False,False
        SB003,S005,.91,False,False,False,False
        SB003,S006,.89,False,False,False,False
        SB004,S007,.89,False,False,False,True
        SB004,S008,.9,False,False,False,True
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
        Group_By_Subject_ID,Sample_ID,Subject_Representative,is_internal_control
        SB001,S001,True,False
        SB001,S002,False,False
        SB002,S003,False,False
        SB002,S004,False,False
        SB002,S005,True,False
        SB003,S006,False,False
        SB003,S007,False,False
        SB004,I001,False,True
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
