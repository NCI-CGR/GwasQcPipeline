import pytest

from cgr_gwas_qc.testing.comparison import file_hashes_equal
from cgr_gwas_qc.workflow.scripts import gtc2plink


################################################################################
# Compare to outputs from old TOP script.
################################################################################
@pytest.mark.regression
def test_gtc2plink_top(fake_data_cache, tmp_path):
    """Run gtc2plink using top strand (default)."""
    gtc2plink.main(
        fake_data_cache / "illumina/gtc/small_genotype.gtc",
        fake_data_cache / "illumina/bpm/small_manifest.bpm",
        "T0001",
        tmp_path,
        "top",
        None,
    )

    expected_map = fake_data_cache / "plink/top/T0001.map"
    expected_ped = fake_data_cache / "plink/top/T0001.ped"

    obs_map = tmp_path / "T0001.map"
    obs_ped = tmp_path / "T0001.ped"

    # Make sure observed and expected results are exactly identical
    assert file_hashes_equal(expected_map, obs_map)
    assert file_hashes_equal(expected_ped, obs_ped)


################################################################################
# Compare to outputs from old FWD script.
################################################################################
@pytest.mark.regression
def test_gtc2plink_fwd(fake_data_cache, tmp_path):
    """Run gtc2plink using fwd strand."""
    gtc2plink.main(
        fake_data_cache / "illumina/gtc/small_genotype.gtc",
        fake_data_cache / "illumina/bpm/small_manifest.bpm",
        "T0001",
        tmp_path,
        "fwd",
        None,
    )

    expected_map = fake_data_cache / "plink/fwd/T0001.map"
    expected_ped = fake_data_cache / "plink/fwd/T0001.ped"

    obs_map = tmp_path / "T0001.map"
    obs_ped = tmp_path / "T0001.ped"

    # Make sure observed and expected results are exactly identical
    assert file_hashes_equal(expected_map, obs_map)
    assert file_hashes_equal(expected_ped, obs_ped)


################################################################################
# Make sure that setting a threshold in creases the number of unknowns (0).
################################################################################
@pytest.mark.regression
def test_gtc2plink_threshold(fake_data_cache, tmp_path):
    """Run gtc2plink using a threshold."""

    # Run without threshold
    no_thresh = tmp_path / "no_thresh"
    no_thresh.mkdir()
    gtc2plink.main(
        fake_data_cache / "illumina/gtc/small_genotype.gtc",
        fake_data_cache / "illumina/bpm/small_manifest.bpm",
        "T0001",
        no_thresh,
        "top",
        None,
    )

    no_thresh_ped = no_thresh / "T0001.ped"

    # Run with threshold
    with_thresh = tmp_path / "with_thresh"
    with_thresh.mkdir()
    gtc2plink.main(
        fake_data_cache / "illumina/gtc/small_genotype.gtc",
        fake_data_cache / "illumina/bpm/small_manifest.bpm",
        "T0001",
        with_thresh,
        "top",
        0.5,
    )

    with_thresh_ped = with_thresh / "T0001.ped"

    # Make sure there are more unknowns (0) when setting a threshold.
    assert no_thresh_ped.read_text().count("0") < with_thresh_ped.read_text().count("0")
