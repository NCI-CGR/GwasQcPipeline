from pathlib import Path

from typer.testing import CliRunner

from cgr_gwas_qc import load_config
from cgr_gwas_qc.cli.pre_flight import _check_user_files, app
from cgr_gwas_qc.testing import chdir

runner = CliRunner()


def test_preflight_refs_only(test_data):
    # GIVEN: a working directory with
    #   - sample sheet
    #   - reference files (BPM, VCF, TBI)
    #   - user files (GTC entry point)
    #   - config

    with chdir(test_data.working_dir):
        # WHEN: Run `cgr pre-flight --no-user-files-check`
        results = runner.invoke(app, ["--user-files-check"])
        # THEN: pre-flight validation passes
        assert results.exit_code == 0


def test_preflight_refs_and_gtc_only(test_data):
    # GIVEN: a working directory with
    #   - sample sheet
    #   - reference files (BPM, VCF, TBI)
    #   - user files (GTC entry point)
    #   - config
    with chdir(test_data.working_dir):
        # WHEN: Run `cgr pre-flight`
        results = runner.invoke(app)
        # THEN: pre-flight validation passes
        assert results.exit_code == 0


def test_check_user_files_no_problems(test_data):
    with chdir(test_data.working_dir):
        cfg = load_config()
        user_files = cfg.config.user_files
        record = dict(cfg.ss.iloc[0, :])
        problems = [problem for problem in _check_user_files(user_files, record) if problem]
        assert 0 == len(problems)


def test_check_user_files_missing_idat_red(test_data):
    with chdir(test_data.working_dir):
        cfg = load_config()
        user_files = cfg.config.user_files
        record = dict(cfg.ss.iloc[0, :])
        Path(user_files.idat_pattern.red.format(**record)).unlink()
        problems = [problem for problem in _check_user_files(user_files, record) if problem]
        assert 1 == len(problems)
        assert "FileNotFound" == problems[0].reason
        assert "idat_red" == problems[0].file_type


def test_check_user_files_missing_gtc(test_data):
    with chdir(test_data.working_dir):
        cfg = load_config()
        user_files = cfg.config.user_files
        record = dict(cfg.ss.iloc[0, :])
        Path(user_files.gtc_pattern.format(**record)).unlink()
        problems = [problem for problem in _check_user_files(user_files, record) if problem]
        assert 1 == len(problems)
        assert "FileNotFound" == problems[0].reason
        assert "gtc" == problems[0].file_type


def test_check_user_files_missing_both_idat(test_data):
    with chdir(test_data.working_dir):
        cfg = load_config()
        user_files = cfg.config.user_files
        record = dict(cfg.ss.iloc[0, :])
        Path(user_files.idat_pattern.red.format(**record)).unlink()
        Path(user_files.idat_pattern.green.format(**record)).unlink()
        problems = [problem for problem in _check_user_files(user_files, record) if problem]
        assert 2 == len(problems)
        assert "FileNotFound" == problems[0].reason
        assert problems[0].file_type in ["idat_red", "idat_green"]
        assert "FileNotFound" == problems[1].reason
        assert problems[1].file_type in ["idat_red", "idat_green"]
