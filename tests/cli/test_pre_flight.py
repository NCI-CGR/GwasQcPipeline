from typer.testing import CliRunner

from cgr_gwas_qc.cli.pre_flight import app
from cgr_gwas_qc.testing import chdir

runner = CliRunner()


def test_preflight_refs_only(gtc_working_dir):
    with chdir(gtc_working_dir):
        results = runner.invoke(app, ["--no-gtc-check", "--no-idat-check"])
        assert results.exit_code == 0


def test_preflight_refs_and_gtc_only(gtc_working_dir):
    with chdir(gtc_working_dir):
        results = runner.invoke(app, ["--gtc-check", "--no-idat-check"])
        assert results.exit_code == 0
