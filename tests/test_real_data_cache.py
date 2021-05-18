"""A number of sanity checks with the real data cache."""
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal, assert_series_equal

from cgr_gwas_qc.parsers.bim import BimFile
from cgr_gwas_qc.testing import file_hashes_equal


################################################################################
# Cached Legacy vs Cached Development regression tests
################################################################################
@pytest.mark.skip(reason="I don't have a BED parser to do the comparisons.")
@pytest.mark.parametrize(
    "legacy_file,dev_file",
    [
        ("legacy_outputs/plink_start/samples.bed", "dev_outputs/sample_level/samples.bed"),
        (
            "legacy_outputs/plink_filter_call_rate_1/samples.bed",
            "dev_outputs/sample_level/call_rate_1/samples.bed",
        ),
        (
            "legacy_outputs/plink_filter_call_rate_2/samples.bed",
            "dev_outputs/sample_level/call_rate_2/samples.bed",
        ),
    ],
)
@pytest.mark.real_data
@pytest.mark.regression
def test_initial_samples_bed(real_data_cache, legacy_file, dev_file):
    """BEDs are not byte perfect matches.

    I would have to do a record-by-record comparisons, but I don't currently
    have a BED parser.
    """
    # legacy_bed = real_data_cache / legacy_file
    # dev_bed = real_data_cache / dev_file
    raise NotImplementedError


@pytest.mark.parametrize(
    "legacy_file,dev_file",
    [
        ("legacy_outputs/plink_start/samples.bim", "dev_outputs/sample_level/samples.bim"),
        (
            "legacy_outputs/plink_filter_call_rate_1/samples.bim",
            "dev_outputs/sample_level/call_rate_1/samples.bim",
        ),
        (
            "legacy_outputs/plink_filter_call_rate_2/samples.bim",
            "dev_outputs/sample_level/call_rate_2/samples.bim",
        ),
    ],
)
@pytest.mark.real_data
@pytest.mark.regression
def test_cached_bims(real_data_cache, legacy_file, dev_file):
    """BIMs are not byte perfect matches.

    The allele order is not guaranteed when generating these files. So I need to
    do a row-by-row comparison.
    """
    legacy_bim = BimFile(real_data_cache / legacy_file)
    dev_bim = BimFile(real_data_cache / dev_file)

    rows_are_equal = [
        legacy == dev for legacy, dev in zip(legacy_bim, dev_bim)
    ]  # NOTE: this ignores allele order
    assert all(rows_are_equal)


@pytest.mark.parametrize(
    "legacy_file,dev_file",
    [
        ("legacy_outputs/plink_start/samples.fam", "dev_outputs/sample_level/samples.fam"),
        (
            "legacy_outputs/plink_filter_call_rate_1/samples.fam",
            "dev_outputs/sample_level/call_rate_1/samples.fam",
        ),
        (
            "legacy_outputs/plink_filter_call_rate_2/samples.fam",
            "dev_outputs/sample_level/call_rate_2/samples.fam",
        ),
    ],
)
@pytest.mark.real_data
@pytest.mark.regression
def test_cached_fams(real_data_cache, legacy_file, dev_file):
    """FAMs should be byte perfect matches."""
    legacy_fam = real_data_cache / legacy_file
    dev_fam = real_data_cache / dev_file
    assert file_hashes_equal(legacy_fam, dev_fam)


@pytest.mark.real_data
@pytest.mark.regression
def test_cached_abf(real_data_cache):
    legacy_abf = real_data_cache / "legacy_outputs/GSAMD-24v1-0_20011747_A1.AF.abf.txt"
    dev_abf = real_data_cache / "dev_outputs/sample_level/GSAMD-24v1-0_20011747_A1.AF.abf.txt"

    legacy = pd.read_csv(legacy_abf, sep="\t")
    dev = pd.read_csv(dev_abf, sep="\t").fillna(0)  # legacy uses 0.0 instead of NA.

    # There were some logic changes so these file are moderately different.
    prop_very_different = (abs(legacy.ABF - dev.ABF) > 1e-6).mean()
    assert 0.01 > prop_very_different  # less than 1% are more than 0.000001 different


@pytest.mark.real_data
@pytest.mark.regression
def test_cached_median_idat_intensity(real_data_cache):
    legacy_file = real_data_cache / "legacy_outputs/all_sample_idat_intensity/idat_intensity.csv"
    dev_file = real_data_cache / "dev_outputs/sample_level/contamination/median_idat_intensity.csv"

    legacy = (
        pd.read_csv(legacy_file)
        .rename(
            {
                "SampId": "Sample_ID",
                "ChipId": "Chip_ID",
                "MedianIntensity": "median_intensity",
            },
            axis=1,
        )
        .set_index("Sample_ID")
    )
    dev = pd.read_csv(dev_file).set_index("Sample_ID")

    assert_frame_equal(legacy, dev, check_like=True)


@pytest.mark.real_data
@pytest.mark.regression
def test_cached_contam_verifiedIDintensity(real_data_cache):
    legacy_file = real_data_cache / "legacy_outputs/all_contam/contam.csv"
    dev_file = real_data_cache / "dev_outputs/sample_level/contamination/verifyIDintensity.csv"

    legacy = (
        pd.read_csv(legacy_file)
        .rename({"ID": "Sample_ID"}, axis=1)
        .set_index("Sample_ID")
        .sort_index()
    )
    dev = pd.read_csv(dev_file).set_index("Sample_ID").sort_index()

    mask = legacy["%Mix"].notna() & dev["%Mix"].notna()

    assert_series_equal(legacy[mask]["%Mix"], dev[mask]["%Mix"], atol=1e-3)
    assert_series_equal(legacy[mask]["LLK"], dev[mask]["LLK"], atol=1000)
    assert_series_equal(legacy[mask]["LLK0"], dev[mask]["LLK0"], atol=500)


################################################################################
# Cached Development vs Current Development regression tests
################################################################################
