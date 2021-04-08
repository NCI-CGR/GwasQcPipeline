"""Test IDAT validation.

IDAT validation is not well developed. It is based on me looking at the
binary output and trying to get a sense of what makes a good IDAT file.
"""
import pytest

from cgr_gwas_qc.testing.data import FakeData
from cgr_gwas_qc.validators.idat import IdatMagicNumberError, IdatTruncatedFileError, validate


@pytest.fixture
def idat_file():
    """Returns the path to a test idat file."""
    return FakeData._data_path / FakeData._test_idat


@pytest.fixture(params=[50, 100])
def truncated_idat_file(tmp_path, idat_file, request):
    """Returns the path to a truncated IDAT file.

    Removes the last ``request.param`` bytes from a test IDAT file and save
    it. Then returns the path to this file.
    """
    data = idat_file.read_bytes()
    trunc_file = tmp_path / "truncated.idat"
    with trunc_file.open("wb") as fh:
        fh.write(data[: -request.param])
    return trunc_file


################################################################################
# Good IDAT file validates
################################################################################
def test_idat_good_idat_file(idat_file):
    # GIVEN: a good IDAT file
    # WHEN-THEN: we validate the file it raises no errors
    validate(idat_file)


################################################################################
# Error if not a IDAT file (IdatMagicNumberError)
################################################################################
def test_IDAT_bad_magic_number():
    # GIVEN: a good bpm file
    # WHEN-THEN: we validate it as a IDAT file it raises an error
    with pytest.raises(IdatMagicNumberError):
        validate(FakeData._data_path / FakeData._illumina_manifest_file)


################################################################################
# Error if file is truncated (IdatTruncatedFileError)
################################################################################
def test_truncated_file(truncated_idat_file):
    # GIVEN: a truncated IDAT file missing a 50+ bytes from the end
    # WHEN-THEN: we validate the file it raises a truncation error
    with pytest.raises(IdatTruncatedFileError):
        validate(truncated_idat_file)
