from pathlib import Path

import pytest

from cgr_gwas_qc.cli import pre_flight
from cgr_gwas_qc.config import ConfigMgr
from cgr_gwas_qc.testing import chdir


def test_preflight(fake_cfg: ConfigMgr):
    # GIVEN: a working directory
    with chdir(fake_cfg.root):
        # WHEN: Run `cgr pre-flight`
        pre_flight.main(Path("config.yml"))
        # THEN: No errors and `cgr_sample_sheet.csv` exists
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
