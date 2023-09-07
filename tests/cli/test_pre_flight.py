from pathlib import Path

import pandas as pd
import pytest
from typer.testing import CliRunner

from cgr_gwas_qc.cli import pre_flight
from cgr_gwas_qc.config import ConfigMgr
from cgr_gwas_qc.testing import chdir

runner = CliRunner()


def test_preflight(fake_cfg: ConfigMgr):
    # GIVEN: a working directory
    with chdir(fake_cfg.root):
        # WHEN: Run `cgr pre-flight`
        res = runner.invoke(pre_flight.app)
        # THEN: No errors and `cgr_sample_sheet.csv` exists
        assert res.exit_code == 0
        assert Path("cgr_sample_sheet.csv").exists()


def test_check_config(fake_cfg: ConfigMgr, capsys):
    # GIVEN: a working directory
    with chdir(fake_cfg.root):
        # WHEN: Run `cgr pre-flight`
        pre_flight.check_config(Path("config.yml"))
        # THEN: No errors and `cgr_sample_sheet.csv` exists
        stdout, _ = capsys.readouterr()
        assert "OK" in stdout


def test_check_config_missing(fake_cfg: ConfigMgr, capsys):
    # GIVEN: a working directory
    with chdir(fake_cfg.root):
        # WHEN: Run `cgr pre-flight`
        with pytest.raises(SystemExit):
            pre_flight.check_config(Path("config_missing.yml"))
        # THEN: No errors and `cgr_sample_sheet.csv` exists
        stdout, _ = capsys.readouterr()
        assert "Config ERROR" in stdout


def test_check_sample_sheet(fake_cfg: ConfigMgr, capsys):
    # GIVEN: a working directory
    with chdir(fake_cfg.root):
        # WHEN: Run `cgr pre-flight`
        pre_flight.check_sample_sheet(
            fake_cfg.config.sample_sheet,
            fake_cfg.config.workflow_params.subject_id_column,
            fake_cfg.config.workflow_params.expected_sex_column,
            fake_cfg.config.workflow_params.case_control_column,
        )
        # THEN: No errors and `cgr_sample_sheet.csv` exists
        stdout, _ = capsys.readouterr()
        assert "OK" in stdout


def test_check_sample_sheet_missing(fake_cfg: ConfigMgr, capsys):
    # GIVEN: a working directory
    with chdir(fake_cfg.root):
        # WHEN: Run `cgr pre-flight`
        with pytest.raises(SystemExit):
            pre_flight.check_sample_sheet(
                Path("sample_sheet_missing.csv"),
                fake_cfg.config.workflow_params.subject_id_column,
                fake_cfg.config.workflow_params.expected_sex_column,
                fake_cfg.config.workflow_params.case_control_column,
            )
        # THEN: No errors and `cgr_sample_sheet.csv` exists
        stdout, _ = capsys.readouterr()
        assert "ERROR" in stdout
        assert "FileNotFound" in stdout


def test_check_reference_files(fake_cfg: ConfigMgr, capsys):
    # GIVEN: a working directory
    with chdir(fake_cfg.root):
        # WHEN: Run `cgr pre-flight`
        pre_flight.check_reference_files(fake_cfg.config.reference_files)
        # THEN: No errors and `cgr_sample_sheet.csv` exists
        stdout, _ = capsys.readouterr()
        assert 3 == stdout.count("OK")


def test_check_reference_files_bpm_error(fake_cfg: ConfigMgr, capsys):
    # GIVEN: a working directory
    with chdir(fake_cfg.root):
        # WHEN: Run `cgr pre-flight`
        reference_files = fake_cfg.config.reference_files.copy()
        reference_files.illumina_manifest_file = "missing.bpm"
        pre_flight.check_reference_files(reference_files)
        # THEN: No errors and `cgr_sample_sheet.csv` exists
        stdout, _ = capsys.readouterr()
        assert "BPM ERROR" in stdout


def test_check_reference_files_vcf_error(fake_cfg: ConfigMgr, capsys):
    # GIVEN: a working directory
    with chdir(fake_cfg.root):
        # WHEN: Run `cgr pre-flight`
        reference_files = fake_cfg.config.reference_files.copy()
        reference_files.thousand_genome_vcf = "missing.vfc"
        pre_flight.check_reference_files(reference_files)
        # THEN: No errors and `cgr_sample_sheet.csv` exists
        stdout, _ = capsys.readouterr()
        assert "VCF ERROR" in stdout


def test_check_reference_files_tbi_error(fake_cfg: ConfigMgr, capsys):
    # GIVEN: a working directory
    with chdir(fake_cfg.root):
        # WHEN: Run `cgr pre-flight`
        reference_files = fake_cfg.config.reference_files.copy()
        reference_files.thousand_genome_tbi = "missing.vfc.tbi"
        pre_flight.check_reference_files(reference_files)
        # THEN: No errors and `cgr_sample_sheet.csv` exists
        stdout, _ = capsys.readouterr()
        assert "VCF.TBI ERROR" in stdout


def test_check_user_files_no_problems(fake_cfg: ConfigMgr):
    with chdir(fake_cfg.root):
        user_files = fake_cfg.config.user_files
        ss = fake_cfg.ss.copy()

        record = dict(ss.iloc[0, :])
        problems = [
            problem for problem in pre_flight._check_user_files(user_files, record) if problem
        ]

        assert 0 == len(problems)


def test_check_user_files_missing_gtc(fake_cfg: ConfigMgr):
    with chdir(fake_cfg.root):
        user_files = fake_cfg.config.user_files
        ss = fake_cfg.ss.copy()

        # Create a record that won't have a GTC
        record = dict(ss.iloc[0, :])
        record["Sample_ID"] = "SB999999"
        problems = [
            problem for problem in pre_flight._check_user_files(user_files, record) if problem
        ]

        assert 1 == len(problems)
        assert "FileNotFound" == problems[0].reason
        assert "gtc" == problems[0].file_type


def test_check_user_files_missing_both_idat(fake_cfg: ConfigMgr):
    with chdir(fake_cfg.root):
        user_files = fake_cfg.config.user_files
        ss = fake_cfg.ss.copy()

        # Create a record that won't have IDAT files
        record = dict(ss.iloc[0, :])
        record["SentrixBarcode_A"] = "999999"
        problems = [
            problem for problem in pre_flight._check_user_files(user_files, record) if problem
        ]

        assert 2 == len(problems)
        assert "FileNotFound" == problems[0].reason
        assert problems[0].file_type in ["idat_red", "idat_green"]
        assert "FileNotFound" == problems[1].reason
        assert problems[1].file_type in ["idat_red", "idat_green"]


################################################################################
# Adding of custom columns based on config options
################################################################################
@pytest.fixture
def updated_sample_sheet(fake_cfg, workflow_params) -> pd.DataFrame:
    ss = fake_cfg.ss.copy()
    return pre_flight.update_sample_sheet(
        ss,
        "LIMS_Individual_ID",
        workflow_params.expected_sex_column,
        workflow_params.case_control_column,
        [
            pre_flight.ProblemFile("SP00005", "test", "test", "test"),
            pre_flight.ProblemFile("SP00003", "FileNotFound", "idat_green", "test"),
        ],
        2,
    )


def test_update_sample_sheet_sample_no_reps(updated_sample_sheet: pd.DataFrame):
    sr = updated_sample_sheet.query("Sample_ID == 'SP00001'").squeeze()
    assert 1 == sr.num_samples_per_subject
    assert pd.isna(sr.replicate_ids)
    assert "M" == sr.expected_sex
    assert "Case" == sr.case_control
    assert not sr.is_internal_control
    assert not sr.is_sample_exclusion


def test_update_sample_sheet_sample_w_reps(updated_sample_sheet: pd.DataFrame):
    sr = updated_sample_sheet.query("Sample_ID == 'SP00002'").squeeze()
    assert 2 == sr.num_samples_per_subject
    assert "SP00002|SP00003" == sr.replicate_ids
    assert "F" == sr.expected_sex
    assert "Control" == sr.case_control
    assert not sr.is_internal_control
    assert not sr.is_sample_exclusion


def test_update_sample_sheet_internal_control(updated_sample_sheet: pd.DataFrame):
    sr = updated_sample_sheet.query("Sample_ID == 'SP00006'").squeeze()
    assert 1 == sr.num_samples_per_subject
    assert pd.isna(sr.replicate_ids)
    assert "M" == sr.expected_sex
    assert "QC" == sr.case_control
    assert sr.is_internal_control
    assert not sr.is_sample_exclusion


def test_update_sample_sheet_excluded_sample(updated_sample_sheet: pd.DataFrame):
    sr = updated_sample_sheet.query("Sample_ID == 'SP00005'").squeeze()
    assert 1 == sr.num_samples_per_subject
    assert pd.isna(sr.replicate_ids)
    assert "M" == sr.expected_sex
    assert "Control" == sr.case_control
    assert not sr.is_internal_control
    assert sr.is_sample_exclusion


def test_update_sample_sheet_missing_idats(updated_sample_sheet: pd.DataFrame):
    assert not updated_sample_sheet.query("Sample_ID == 'SP00005'").is_missing_idats.squeeze()
    assert updated_sample_sheet.query("Sample_ID == 'SP00003'").is_missing_idats.squeeze()
