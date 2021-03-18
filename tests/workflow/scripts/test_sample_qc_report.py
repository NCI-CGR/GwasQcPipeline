from io import StringIO
from pathlib import Path
from textwrap import dedent

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_series_equal

from cgr_gwas_qc import load_config
from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.testing.data import RealData
from cgr_gwas_qc.workflow.scripts.sample_qc_table import CASE_CONTROL_DTYPE, SEX_DTYPE


@pytest.mark.real_data
@pytest.fixture(scope="session")
def data_cache_and_cfg(tmp_path_factory, pytestconfig):
    """Real data and config object for QC table unit tests."""
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    session_tmp_path: Path = tmp_path_factory.mktemp("sample_qc_table_script")
    data_cache = (
        RealData(session_tmp_path, full_sample_sheet=False)
        .copy_sample_sheet()
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
    )

    with chdir(session_tmp_path):
        cfg = load_config()

    return data_cache, cfg


# @pytest.mark.skip(reason="Flaky test, need to look into it")
@pytest.mark.real_data
@pytest.mark.parametrize("expected_sex_col", ["Expected_Sex", "Identifiler_Sex"])
def test_wrangle_sample_sheet(data_cache_and_cfg, expected_sex_col):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _wrangle_sample_sheet

    # GIVEN: A test sample sheet and column to use as `expected_sex_col`
    _, cfg = data_cache_and_cfg

    # WHEN: I wrangle the sample sheet
    ss = _wrangle_sample_sheet(cfg.ss, expected_sex_col)

    # THEN: Basic properties
    assert ss.index.name == "Sample_ID"
    assert isinstance(ss.internal_control.dtype, pd.BooleanDtype)
    assert ss.expected_sex.dtype == SEX_DTYPE
    assert ss.case_control.dtype == CASE_CONTROL_DTYPE

    # The `expected_sex` column should be the same as the column passed as `expected_sex_col`
    assert_series_equal(
        ss.expected_sex.sort_index().astype("object"),
        cfg.ss.set_index("Sample_ID")[expected_sex_col].sort_index(),
        check_names=False,
    )


@pytest.mark.real_data
def test_read_imiss_start(data_cache_and_cfg):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_imiss_start

    # GIVEN: A test sample sheet, real production outputs, and a list of Sample_IDs
    data_cache, cfg = data_cache_and_cfg
    file_name = data_cache / "production_outputs/plink_start/samples_start.imiss"
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: I parse the imiss table.
    sr = _read_imiss_start(file_name, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "Call_Rate_Initial"


@pytest.mark.real_data
def test_read_imiss_cr1(data_cache_and_cfg):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_imiss_cr1

    # GIVEN: A test sample sheet, real production outputs, and a list of Sample_IDs
    data_cache, cfg = data_cache_and_cfg
    file_name = data_cache / "production_outputs/plink_filter_call_rate_1/samples_filter1.imiss"
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: I parse the imiss table.
    df = _read_imiss_cr1(file_name, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "Sample_ID"
    assert "Call_Rate_1" in df.columns
    assert "Call_Rate_1_filter" in df.columns


@pytest.mark.real_data
def test_read_imiss_cr2(data_cache_and_cfg):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_imiss_cr2

    # GIVEN: A test sample sheet, real production outputs, and a list of Sample_IDs
    data_cache, cfg = data_cache_and_cfg
    file_name = data_cache / "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss"
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: I parse the imiss table.
    df = _read_imiss_cr2(file_name, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "Sample_ID"
    assert "Call_Rate_2" in df.columns
    assert "Call_Rate_2_filter" in df.columns


@pytest.mark.real_data
def test_read_sexcheck_cr1(data_cache_and_cfg):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_sexcheck_cr1

    # GIVEN: A test sample sheet, real production outputs, and the expected sex
    # calls from the sample table
    data_cache, cfg = data_cache_and_cfg
    file_name = data_cache / "production_outputs/plink_filter_call_rate_1/samples_filter1.sexcheck"
    expected_sex_calls = cfg.ss.set_index("Sample_ID")["Expected_Sex"]
    expected_sex_calls[
        "SB_missing_sample"
    ] = "F"  # add an extra sample to check how missing values are handled

    # WHEN: read the sexcheck table
    df = _read_sexcheck_cr1(file_name, expected_sex_calls)

    # THEN: Basic properties
    assert df.index.name == "Sample_ID"
    assert df.Predicted_Sex.dtype == pd.CategoricalDtype(categories=["M", "F", "U"])
    assert isinstance(df.sex_discordant.dtype, pd.BooleanDtype)
    assert df.ChrX_Inbreed_estimate.dtype is np.dtype("float64")


@pytest.mark.real_data
def test_read_ancestry_GRAF(data_cache_and_cfg):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_ancestry

    # GIVEN: A test sample sheet, example outputs from GRAF -pop, and a list of Sample_IDs
    _, cfg = data_cache_and_cfg
    ancestry = cfg.root / "graf.tsv"
    ancestry.write_text(
        dedent(
            """
            Subject\tSelf-reported ancestry\tGD1\tGD2\tGD3\tGD4\tP_f (%)\tP_e (%)\tP_a (%)\tPopID\tComputed population
            SB034307_PC26469_B09\t\t1.63\t1.45\t0.81\t0.03\t0.6\t99.1\t0.2\t1\tEuropean
            SB034327_PC26469_C09\t\t1.64\t1.47\t0.85\t0.02\t0.0\t100.0\t0.0\t1\tEuropean
            """
        )
    )
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: Parse the GRAF output
    df = _read_ancestry(ancestry, Sample_IDs)

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
def test_read_known_replicates(data_cache_and_cfg):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_known_replicates

    # GIVEN: A test sample sheet, config, real production outputs, and a list of Sample_IDs
    data_cache, cfg = data_cache_and_cfg
    cutoff = cfg.config.software_params.dup_concordance_cutoff
    reps = data_cache / "production_outputs/concordance/KnownReplicates.csv"
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: Check for discordant replicates
    sr = _read_known_replicates(reps, cutoff, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "Expected Replicate Discordance"
    assert sr.dtype is np.dtype("bool")


@pytest.mark.real_data
def test_read_unknown_replicates(data_cache_and_cfg):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_unknown_replicates

    # GIVEN: A test sample sheet, real production outputs, and a list of Sample_IDs
    data_cache, cfg = data_cache_and_cfg
    reps = data_cache / "production_outputs/concordance/UnknownReplicates.csv"
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: I parse unknown concordant samples table.
    sr = _read_unknown_replicates(reps, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "Unexpected Replicate"
    assert sr.dtype is np.dtype("bool")


@pytest.mark.real_data
def test_read_contam(data_cache_and_cfg):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_contam

    # GIVEN: A test sample sheet, config, real production outputs, and a list of Sample_IDs
    data_cache, cfg = data_cache_and_cfg
    cutoff = cfg.config.software_params.contam_threshold
    contam = data_cache / "production_outputs/all_contam/contam.csv"
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: Parse the contamination table
    df = _read_contam(contam, cutoff, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "Sample_ID"
    assert "Contamination_Rate" in df.columns
    assert "Contaminated" in df.columns
    assert df.dtypes["Contaminated"] is np.dtype("bool")
    assert df.dtypes["Contamination_Rate"] is np.dtype("float64")


@pytest.mark.real_data
def test_read_contam_file_name_none(data_cache_and_cfg):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_contam

    # GIVEN: A test sample sheet, config, real production outputs, and a list of Sample_IDs
    _, cfg = data_cache_and_cfg
    cutoff = cfg.config.software_params.contam_threshold
    contam = None
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: Parse the contamination table.
    df = _read_contam(contam, cutoff, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "Sample_ID"
    assert "Contamination_Rate" in df.columns
    assert "Contaminated" in df.columns


@pytest.mark.real_data
def test_read_intensity(data_cache_and_cfg):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_intensity

    # GIVEN: A test sample sheet, real production outputs, and a list of Sample_IDs
    data_cache, cfg = data_cache_and_cfg
    intensity = data_cache / "production_outputs/all_sample_idat_intensity/idat_intensity.csv"
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: Parse the read intensity table.
    sr = _read_intensity(intensity, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "IdatIntensity"
    assert sr.dtype is np.dtype("float64")


@pytest.mark.real_data
def test_read_intensity_file_name_none(data_cache_and_cfg):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _read_intensity

    # GIVEN: A test sample sheet, No production outputs, and a list of Sample_IDs
    _, cfg = data_cache_and_cfg
    intensity = None
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: Parse the read intensity table with file name of None
    sr = _read_intensity(intensity, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "IdatIntensity"


@pytest.mark.real_data
def test_check_idat_files(data_cache_and_cfg):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _check_idats_files

    # GIVEN: A test sample sheet and real production outputs
    _, cfg = data_cache_and_cfg

    # WHEN: Scan the folder for Idat files.
    sr = _check_idats_files(cfg)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "IdatsInProjectDir"

    # All idat files should be found
    assert all(sr)


@pytest.mark.real_data
def test_check_idat_files_one_missing(data_cache_and_cfg):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _check_idats_files

    # GIVEN: A test sample sheet, real production outputs, and a fake sample
    # with missing Idat files.
    _, cfg = data_cache_and_cfg

    # fake entry
    fake = cfg.ss.iloc[0, :].copy()
    fake["Sample_ID"] = "fake_Sample_ID"
    fake["SentrixBarcode_A"] = "fake_barcode"
    cfg.ss = cfg.ss.append(fake, ignore_index=True)

    # WHEN: Scan the folder for Idat files.
    sr = _check_idats_files(cfg)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "IdatsInProjectDir"

    # I should have 2 samples with idat files found and 1 missing these files.
    assert sum(sr) == 2
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
        Group_By_Subject_ID,Sample_ID,Call_Rate_2,Low Call Rate,Contaminated,Expected Replicate Discordance,internal_control
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

    # WHEN: I loook for the representative sample for each Subject Group
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
        Group_By_Subject_ID,Sample_ID,Subject_Representative,internal_control
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
