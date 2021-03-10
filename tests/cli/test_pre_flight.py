from typer.testing import CliRunner

from cgr_gwas_qc.cli.pre_flight import app
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
