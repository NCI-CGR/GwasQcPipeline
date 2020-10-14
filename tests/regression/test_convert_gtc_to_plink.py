from pathlib import Path

import pytest
from typer.testing import CliRunner

from cgr_gwas_qc.testing import file_hashes_equal
from cgr_gwas_qc.workflow.scripts.gtc2plink import app

runner = CliRunner()

sample_id = "T0001"


################################################################################
# Compare to outputs from old TOP script.
################################################################################
@pytest.mark.regression
def test_gtc2plink_top(bpm_file, gtc_file, tmpdir):
    """Run gtc2plink using top strand (default)."""

    exp_map = Path("tests/data/plink/top/T0001.map")
    exp_ped = Path("tests/data/plink/top/T0001.ped")

    obs_map = Path(tmpdir) / f"{sample_id}.map"
    obs_ped = Path(tmpdir) / f"{sample_id}.ped"

    # Run script
    results = runner.invoke(app, [gtc_file.as_posix(), bpm_file.as_posix(), sample_id, str(tmpdir)])
    assert results.exit_code == 0

    # Check STDOUT and make sure it is working
    assert f"Saving MAP file: {obs_map}" in results.stdout
    assert f"Saving PED file: {obs_ped}" in results.stdout

    # Make sure observed and expecte results are exactly identical
    assert file_hashes_equal(obs_map, exp_map)
    assert file_hashes_equal(obs_ped, exp_ped)


################################################################################
# Compare to outputs from old FWD script.
################################################################################
@pytest.mark.regression
def test_gtc2plink_fwd(bpm_file, gtc_file, tmpdir):
    """Run gtc2plink using fwd strand."""
    exp_map = Path("tests/data/plink/fwd/T0001.map")
    exp_ped = Path("tests/data/plink/fwd/T0001.ped")

    obs_map = Path(tmpdir) / f"{sample_id}.map"
    obs_ped = Path(tmpdir) / f"{sample_id}.ped"

    # Run script
    results = runner.invoke(
        app, [gtc_file.as_posix(), bpm_file.as_posix(), sample_id, str(tmpdir), "--strand", "FWD"]
    )
    assert results.exit_code == 0

    # Check STDOUT and make sure it is working
    assert f"Saving MAP file: {obs_map}" in results.stdout
    assert f"Saving PED file: {obs_ped}" in results.stdout

    # Make sure observed and expecte results are exactly identical
    assert file_hashes_equal(obs_map, exp_map)
    assert file_hashes_equal(obs_ped, exp_ped)


################################################################################
# Make sure that setting a threshold in creases the number of unknowns (0).
################################################################################
@pytest.mark.regression
def test_gtc2plink_threshold(bpm_file, gtc_file, tmpdir):
    """Run gtc2plink using a genotype call threshold."""

    # Run without threshold
    ped_no_thresh = Path(tmpdir) / f"{sample_id}.ped"

    results = runner.invoke(app, [gtc_file.as_posix(), bpm_file.as_posix(), sample_id, str(tmpdir)])
    assert results.exit_code == 0

    # Run with threshold
    threshold = 0.5
    _sample_id = f"{sample_id}_{threshold}"
    ped_thresh = Path(tmpdir) / f"{_sample_id}.ped"

    results = runner.invoke(
        app,
        [gtc_file.as_posix(), bpm_file.as_posix(), _sample_id, str(tmpdir), "--cutoff", threshold],
    )
    assert results.exit_code == 0

    # Make sure there are more unknowns (0) when setting a threshold.
    assert ped_no_thresh.read_text().count("0") < ped_thresh.read_text().count("0")
