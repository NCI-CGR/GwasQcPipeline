from typer.testing import CliRunner

from cgr_gwas_qc.cli.pre_flight import app
from cgr_gwas_qc.testing import chdir

runner = CliRunner()


def test_preflight_refs_only(working_dir):
    with chdir(working_dir):
        results = runner.invoke(app, ["--no-gtc-check", "--no-idat-check"])
        assert results.exit_code == 0


def test_preflight_refs_and_gtc_only(working_dir):
    with chdir(working_dir):
        results = runner.invoke(app, ["--gtc-check", "--no-idat-check"])
        assert results.exit_code == 0
