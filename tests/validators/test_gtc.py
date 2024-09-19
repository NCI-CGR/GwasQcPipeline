"""Test GTC validation.

Validation is provided by the Illumina parsing library. Since this is a
binary format, it would be really hard to test individual issues. For now I
am just testing that a good file passes.
"""

import pytest

from cgr_gwas_qc.testing.data import FakeData
from cgr_gwas_qc.validators.gtc import GtcMagicNumberError, GtcTruncatedFileError, validate


@pytest.fixture
def gtc_file():
    """Returns the path to a test idat file."""
    return FakeData._data_path / FakeData._test_gtc


@pytest.fixture(params=[1, 2, 4])
def truncated_gtc_file(tmp_path, gtc_file, request):
    """Returns the path to a truncated GTC file.

    Removes the last ``request.param`` bytes from a test GTC file and save
    it. Then returns the path to this file.
    """
    data = gtc_file.read_bytes()
    trunc_file = tmp_path / "truncated.gtc"
    with trunc_file.open("wb") as fh:
        fh.write(data[: -request.param])
    return trunc_file


################################################################################
# Good GTC file validates
################################################################################
def test_gtc_good_gtc_file(gtc_file):
    # GIVEN: a good gtc file
    # WHEN-THEN: we validate the file it raises no errors
    validate(gtc_file)


################################################################################
# Error if not a GTC file (GtcMagicNumberError)
################################################################################
def test_gtc_bad_magic_number():
    # GIVEN: a good bpm file
    # WHEN-THEN: we validate it as a GTC file it raises an error
    with pytest.raises(GtcMagicNumberError):
        validate(FakeData._data_path / FakeData._illumina_manifest_file)


################################################################################
# Error if file is truncated (GtcTruncatedFileError)
################################################################################
def test_truncated_file(truncated_gtc_file):
    # GIVEN: a truncated GTC file missing a few bytes from the end
    # WHEN-THEN: we validate the file it raises a truncation error
    with pytest.raises(GtcTruncatedFileError):
        validate(truncated_gtc_file)
