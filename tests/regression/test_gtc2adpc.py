from pathlib import Path

import pytest
from numpy import isclose
from typer.testing import CliRunner

from cgr_gwas_qc.cli.gtc2adpc import app
from cgr_gwas_qc.parsers.adpc import AdpcReader

runner = CliRunner()


@pytest.fixture
def expected_adpc():
    return Path("tests/data/small.adpc.bin")


@pytest.fixture
def expected_snp_count():
    return Path("tests/data/small.adpc.bin.numSnps.txt")


@pytest.mark.regression
def test_gtc2adpc(gtc_file, bpm_file, expected_adpc, expected_snp_count, tmpdir):
    obs_adpc = Path(tmpdir) / "adpc.bin"
    obs_snp_count = Path(tmpdir) / "snp_counts.txt"

    results = runner.invoke(
        app,
        [gtc_file.as_posix(), bpm_file.as_posix(), obs_adpc.as_posix(), obs_snp_count.as_posix()],
    )
    assert results.exit_code == 0
    assert expected_snp_count.read_text() == obs_snp_count.read_text()

    with AdpcReader(expected_adpc) as expect, AdpcReader(obs_adpc) as obs:
        for exp_row, obs_row in zip(expect, obs):
            assert exp_row.x_raw == obs_row.x_raw
            assert exp_row.y_raw == obs_row.y_raw
            assert isclose(exp_row.x_norm, obs_row.x_norm)
            assert isclose(exp_row.y_norm, obs_row.y_norm)
            assert isclose(exp_row.genotype_score, obs_row.genotype_score)
            assert exp_row.genotype == obs_row.genotype
