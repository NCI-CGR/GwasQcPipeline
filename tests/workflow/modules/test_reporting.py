from pathlib import Path
from textwrap import dedent
from typing import Tuple

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal, assert_series_equal
from typer.testing import CliRunner

from cgr_gwas_qc import load_config
from cgr_gwas_qc.config import ConfigMgr
from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.testing.data import RealData

runner = CliRunner()


@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_data_cache(tmp_path_factory):
    """Real data and config object for QC table unit tests."""
    session_tmp_path: Path = tmp_path_factory.mktemp("data_for_sample_qc_report")
    data_path = (
        RealData(session_tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
    )

    with chdir(session_tmp_path):
        cfg = load_config()

    return data_path, cfg


@pytest.mark.parametrize("expected_sex_col", ["Expected_Sex", "LIMS_Individual_ID"])
def test_wrangle_sample_sheet(real_data_cache: Tuple[RealData, ConfigMgr], expected_sex_col):
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import _wrangle_sample_sheet

    # GIVEN: A test sample sheet and column to use as `Expected_Sex`
    _, cfg = real_data_cache

    # WHEN: I wrangle the sample sheet
    ss = _wrangle_sample_sheet(cfg.ss, expected_sex_col)

    # THEN: Basic properties
    assert isinstance(ss, pd.DataFrame)
    assert ss.index.name == "Sample_ID"

    # The `Expected_Sex` column should be the same as the column passed as `expected_sex_col`
    assert_series_equal(
        ss["Expected_Sex"], cfg.ss.set_index("Sample_ID")[expected_sex_col], check_names=False
    )


@pytest.mark.real_data
def test_read_imiss_start(real_data_cache: Tuple[RealData, ConfigMgr]):
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import _read_imiss_start

    # GIVEN: A test sample sheet, real production outputs, and a list of Sample_IDs
    data_path, cfg = real_data_cache
    file_name = data_path / "production_outputs/plink_start/samples_start.imiss"
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: I parse the imiss table.
    sr = _read_imiss_start(file_name, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "Call_Rate_Initial"


@pytest.mark.real_data
def test_read_imiss_cr1(real_data_cache: Tuple[RealData, ConfigMgr]):
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import _read_imiss_cr1

    # GIVEN: A test sample sheet, real production outputs, and a list of Sample_IDs
    data_path, cfg = real_data_cache
    file_name = data_path / "production_outputs/plink_filter_call_rate_1/samples_filter1.imiss"
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: I parse the imiss table.
    df = _read_imiss_cr1(file_name, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "Sample_ID"
    assert "Call_Rate_1" in df.columns
    assert "Call_Rate_1_filter" in df.columns


@pytest.mark.real_data
def test_read_imiss_cr2(real_data_cache: Tuple[RealData, ConfigMgr]):
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import _read_imiss_cr2

    # GIVEN: A test sample sheet, real production outputs, and a list of Sample_IDs
    data_path, cfg = real_data_cache
    file_name = data_path / "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss"
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: I parse the imiss table.
    df = _read_imiss_cr2(file_name, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "Sample_ID"
    assert "Call_Rate_2" in df.columns
    assert "Call_Rate_2_filter" in df.columns


@pytest.mark.real_data
def test_read_sexcheck_cr1(real_data_cache: Tuple[RealData, ConfigMgr]):
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import _read_sexcheck_cr1

    # GIVEN: A test sample sheet, real production outputs, and the expected sex
    # calls from the sample table
    data_path, cfg = real_data_cache
    file_name = data_path / "production_outputs/plink_filter_call_rate_1/samples_filter1.sexcheck"
    expected_sex_calls = cfg.ss.set_index("Sample_ID")["Expected_Sex"]

    # WHEN: read the sexcheck table
    df = _read_sexcheck_cr1(file_name, expected_sex_calls)

    # THEN: Basic properties
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "Sample_ID"
    assert "ChrX_Inbreed_estimate" in df.columns
    assert "Predicted_Sex" in df.columns
    assert "SexMatch" in df.columns
    assert "Sex Discordant" in df.columns
    assert df.dtypes["Sex Discordant"] is np.dtype("bool")
    assert df.dtypes["ChrX_Inbreed_estimate"] is np.dtype("float64")


@pytest.mark.real_data
def test_read_ancestry(real_data_cache: Tuple[RealData, ConfigMgr]):
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import _read_ancestry

    # GIVEN: A test sample sheet, example outputs from GRAF -pop, and a list of Sample_IDs
    _, cfg = real_data_cache
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
def test_read_known_replicates(real_data_cache: Tuple[RealData, ConfigMgr]):
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import _read_known_replicates

    # GIVEN: A test sample sheet, config, real production outputs, and a list of Sample_IDs
    data_path, cfg = real_data_cache
    cutoff = cfg.config.software_params.dup_concordance_cutoff
    reps = data_path / "production_outputs/concordance/KnownReplicates.csv"
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: Check for discordant replicates
    sr = _read_known_replicates(reps, cutoff, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "Expected Replicate Discordance"
    assert sr.dtype is np.dtype("bool")


@pytest.mark.real_data
def test_read_unknown_replicates(real_data_cache: Tuple[RealData, ConfigMgr]):
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import _read_unknown_replicates

    # GIVEN: A test sample sheet, real production outputs, and a list of Sample_IDs
    data_path, cfg = real_data_cache
    reps = data_path / "production_outputs/concordance/UnknownReplicates.csv"
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: I parse unknown concordant samples table.
    sr = _read_unknown_replicates(reps, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "Unexpected Replicate"
    assert sr.dtype is np.dtype("bool")


@pytest.mark.real_data
def test_read_contam(real_data_cache: Tuple[RealData, ConfigMgr]):
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import _read_contam

    # GIVEN: A test sample sheet, config, real production outputs, and a list of Sample_IDs
    data_path, cfg = real_data_cache
    cutoff = cfg.config.software_params.contam_threshold
    contam = data_path / "production_outputs/all_contam/contam.csv"
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
def test_read_contam_file_name_none(real_data_cache: Tuple[RealData, ConfigMgr]):
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import _read_contam

    # GIVEN: A test sample sheet, config, real production outputs, and a list of Sample_IDs
    _, cfg = real_data_cache
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
def test_read_intensity(real_data_cache: Tuple[RealData, ConfigMgr]):
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import _read_intensity

    # GIVEN: A test sample sheet, real production outputs, and a list of Sample_IDs
    data_path, cfg = real_data_cache
    intensity = data_path / "production_outputs/all_sample_idat_intensity/idat_intensity.csv"
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: Parse the read intensity table.
    sr = _read_intensity(intensity, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "IdatIntensity"
    assert sr.dtype is np.dtype("float64")


@pytest.mark.real_data
def test_read_intensity_file_name_none(real_data_cache: Tuple[RealData, ConfigMgr]):
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import _read_intensity

    # GIVEN: A test sample sheet, No production outputs, and a list of Sample_IDs
    _, cfg = real_data_cache
    intensity = None
    Sample_IDs = cfg.ss.set_index("Sample_ID").index

    # WHEN: Parse the read intensity table with file name of None
    sr = _read_intensity(intensity, Sample_IDs)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "IdatIntensity"


@pytest.mark.real_data
def test_check_idat_files(real_data_cache: Tuple[RealData, ConfigMgr]):
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import _check_idats_files

    # GIVEN: A test sample sheet and real production outputs
    _, cfg = real_data_cache

    # WHEN: Scan the folder for Idat files.
    sr = _check_idats_files(cfg)

    # THEN: Basic properties
    assert isinstance(sr, pd.Series)
    assert sr.index.name == "Sample_ID"
    assert sr.name == "IdatsInProjectDir"

    # All idat files should be found
    assert all(sr == "YES")


@pytest.mark.real_data
def test_check_idat_files_one_missing(real_data_cache: Tuple[RealData, ConfigMgr]):
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import _check_idats_files

    # GIVEN: A test sample sheet, real production outputs, and a fake sample
    # with missing Idat files.
    _, cfg = real_data_cache

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
    assert sum(sr == "YES") == 2
    assert sum(sr == "NO") == 1


@pytest.mark.regression
@pytest.mark.real_data
def test_sample_qc_report(tmp_path: Path, monkeypatch):
    """Integration test for the Internal Sample QC Report"""
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import app

    # GIVEN: A real data store with the full sample table and config.
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
    )

    # I don't have GRAF output from the production workflow. I am hacking my
    # way around this by:
    # passing an empty file to the script
    graf = tmp_path / "ancestry/graf_ancestry_calls.txt"
    graf.parent.mkdir(parents=True, exist_ok=True)
    graf.touch()

    # monkeypatching the `_read_ancestry` function to just return None which is
    # ignored by `pd.concat`
    def mock_ancestry(*args, **kwargs):
        return None

    monkeypatch.setattr(
        "cgr_gwas_qc.workflow.scripts.sample_qc_report._read_ancestry", mock_ancestry
    )

    # WHEN: I change directory (b/c I use `load_config` in the script) and run the script.
    with chdir(tmp_path):
        results = runner.invoke(
            app,
            [
                (data_store / "production_outputs/plink_start/samples_start.imiss").as_posix(),
                (
                    data_store / "production_outputs/plink_filter_call_rate_1/samples_filter1.imiss"
                ).as_posix(),
                (
                    data_store / "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss"
                ).as_posix(),
                (
                    data_store
                    / "production_outputs/plink_filter_call_rate_1/samples_filter1.sexcheck"
                ).as_posix(),
                graf.as_posix(),  # This is an empty file
                (data_store / "production_outputs/concordance/KnownReplicates.csv").as_posix(),
                (data_store / "production_outputs/concordance/UnknownReplicates.csv").as_posix(),
                "--contam",
                (data_store / "production_outputs/all_contam/contam.csv").as_posix(),
                "--intensity",
                (
                    data_store / "production_outputs/all_sample_idat_intensity/idat_intensity.csv"
                ).as_posix(),
                (tmp_path / "all_sample_qc.csv").as_posix(),
                (tmp_path / "qc_table_for_LimsUpload").as_posix(),
            ],
            catch_exceptions=False,
        )

    # THEN: the script runs successfully
    assert results.exit_code == 0  # Make sure it ran successfully

    # Outputs (excluding Ancestry) should match with legacy workflow. GRAF is
    # not part of the legacy workflow so I am ignoring Ancestry results.
    exclude_cols = ["AFR", "ASN", "EUR", "Ancestry"]  # Ignore GRAF columns
    exclude_cols.append(
        "IdatsInProjectDir"
    )  # Observed values will be wrong b/c Idats were not copied locally
    obs_ = (
        pd.read_csv(tmp_path / "all_sample_qc.csv")
        .drop(exclude_cols, axis=1)
        .sort_values("Sample_ID")
        .reset_index(drop=True)
    )
    exp_ = (
        pd.read_csv(data_store / "production_outputs/all_sample_qc.csv")
        .drop(exclude_cols, axis=1)
        .sort_values("Sample_ID")
        .reset_index(drop=True)
    )
    assert_frame_equal(obs_, exp_)


@pytest.mark.real_data
def test_summary_stats(real_data_cache: Tuple[RealData, ConfigMgr]):
    """Test the createion of summary_stats.txt

    I am not able to perform regression testing because there are some issues
    with the original script. First, they have the contamination threshold of
    0.1 hardcoded so the summary stats don't match the sample_qc_table if
    this threshold is changed by the user. Second, the fail call rate
    includes NAs, where in the sample_qc_table these are converted to True.
    Third, because I switched from R to python, there are some small
    naming/formatting differences that I felt were not necessary to repeat.
    """
    from cgr_gwas_qc.workflow.scripts.sample_qc_report_summary_stats import app

    # GIVEN: A config and the all_sample_qc.csv file.
    data_path, cfg = real_data_cache

    # WHEN: I run in the same directory as the config.yml file.
    with chdir(cfg.root):
        results = runner.invoke(
            app,
            [(data_path / "production_outputs/all_sample_qc.csv").as_posix(), "summary_stats.txt"],
        )

    # THEN: The script runs to completion
    assert results.exit_code == 0
    # The file is created
    assert (cfg.root / "summary_stats.txt").exists()
    # The file is not empty
    assert len((cfg.root / "summary_stats.txt").read_text().splitlines()) == 60
