from pathlib import Path

import pytest
from typer.testing import CliRunner

from cgr_gwas_qc.cli.gtc2plink import app
from cgr_gwas_qc.testing import file_hashes_equal

runner = CliRunner()

sample_id = "T0001"


@pytest.mark.regression
def test_gtc2plink_top(bpm_file, gtc_file, tmpdir):
    """Run gtc2plink using top strand (default)."""
    results = runner.invoke(app, [gtc_file.as_posix(), bpm_file.as_posix(), sample_id, str(tmpdir)])
    assert results.exit_code == 0
    assert f"Saving MAP file: {tmpdir}/{sample_id}.map" in results.stdout
    assert f"Saving PED file: {tmpdir}/{sample_id}.ped" in results.stdout

    obs_map = Path(tmpdir) / f"{sample_id}.map"
    exp_map = Path("tests/data/plink/top/T0001.map")
    assert file_hashes_equal(obs_map, exp_map)

    obs_ped = Path(tmpdir) / f"{sample_id}.ped"
    exp_ped = Path("tests/data/plink/top/T0001.ped")
    assert file_hashes_equal(obs_ped, exp_ped)


@pytest.mark.regression
def test_gtc2plink_fwd(bpm_file, gtc_file, tmpdir):
    """Run gtc2plink using fwd strand."""
    results = runner.invoke(
        app, [gtc_file.as_posix(), bpm_file.as_posix(), sample_id, str(tmpdir), "--strand", "FWD"]
    )
    assert results.exit_code == 0
    assert f"Saving MAP file: {tmpdir}/{sample_id}.map" in results.stdout
    assert f"Saving PED file: {tmpdir}/{sample_id}.ped" in results.stdout

    obs_map = Path(tmpdir) / f"{sample_id}.map"
    exp_map = Path("tests/data/plink/fwd/T0001.map")
    assert file_hashes_equal(obs_map, exp_map)

    obs_ped = Path(tmpdir) / f"{sample_id}.ped"
    exp_ped = Path("tests/data/plink/fwd/T0001.ped")
    assert file_hashes_equal(obs_ped, exp_ped)


@pytest.mark.regression
def test_gtc2plink_threshold(bpm_file, gtc_file, tmpdir):
    """Run gtc2plink using a genotype call threshold."""
    # Run without threshold
    results = runner.invoke(app, [gtc_file.as_posix(), bpm_file.as_posix(), sample_id, str(tmpdir)])
    assert results.exit_code == 0
    ped_no_thresh = Path(tmpdir) / f"{sample_id}.ped"

    # Run with threshold
    threshold = 0.5
    _sample_id = f"{sample_id}_{threshold}"
    results = runner.invoke(
        app,
        [gtc_file.as_posix(), bpm_file.as_posix(), _sample_id, str(tmpdir), "--cutoff", threshold],
    )
    assert results.exit_code == 0
    ped_thresh = Path(tmpdir) / f"{_sample_id}.ped"

    # Then there should be more unknowns (0).
    assert ped_no_thresh.read_text().count("0") < ped_thresh.read_text().count("0")
