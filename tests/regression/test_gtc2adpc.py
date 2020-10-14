from pathlib import Path

import pytest
from numpy import isclose
from typer.testing import CliRunner

from cgr_gwas_qc.parsers.illumina import AdpcReader
from cgr_gwas_qc.workflow.scripts.gtc2adpc import app

runner = CliRunner()


################################################################################
# Compare to outputs from old script.
################################################################################
@pytest.mark.regression
def test_gtc2adpc(gtc_file, bpm_file, tmpdir):
    expected_adpc = Path("tests/data/small.adpc.bin")
    expected_snp_count = Path("tests/data/small.adpc.bin.numSnps.txt")

    obs_adpc = Path(tmpdir) / "adpc.bin"
    obs_snp_count = Path(tmpdir) / "snp_counts.txt"

    # Run the script
    results = runner.invoke(
        app,
        [gtc_file.as_posix(), bpm_file.as_posix(), obs_adpc.as_posix(), obs_snp_count.as_posix()],
    )
    assert results.exit_code == 0

    # Make sure snp counts are identical
    assert expected_snp_count.read_text() == obs_snp_count.read_text()

    # Make sure outputs are about equal. The output contains floats, due to
    # rounding differences between py2 and py3 it is impossible to get the
    # exact same number. However, they should be extremely close.
    with AdpcReader(expected_adpc) as expect, AdpcReader(obs_adpc) as obs:
        for exp_row, obs_row in zip(expect, obs):
            # ints should be the exact same
            assert exp_row.x_raw == obs_row.x_raw
            assert exp_row.y_raw == obs_row.y_raw
            assert exp_row.genotype == obs_row.genotype

            # floats should be really close
            assert isclose(exp_row.x_norm, obs_row.x_norm)
            assert isclose(exp_row.y_norm, obs_row.y_norm)
            assert isclose(exp_row.genotype_score, obs_row.genotype_score)
