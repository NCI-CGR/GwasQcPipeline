"""Test BGZIP validation.

Since this is a binary format, it would be really hard to test individual
issues. For now I am just testing that a good file passes.
"""
import pytest

from cgr_gwas_qc.testing.data import FakeData
from cgr_gwas_qc.validators.bgzip import BgzipMagicNumberError, BgzipTruncatedFileError, validate


@pytest.fixture
def vcf_file():
    """Returns the path to a test idat file."""
    return FakeData._data_path / FakeData._thousand_genome_vcf


@pytest.fixture(params=[1, 2, 4])
def truncated_vcf_file(tmp_path, vcf_file, request):
    """Returns the path to a truncated VCF file.

    Removes the last ``request.param`` bytes from a test VCF file and save
    it. Then returns the path to this file.
    """
    data = vcf_file.read_bytes()
    trunc_file = tmp_path / "truncated.vcf.gz"
    with trunc_file.open("wb") as fh:
        fh.write(data[: -request.param])
    return trunc_file


################################################################################
# Good VCF and TBI files validate ok
################################################################################
def test_bgzip_good_vcf_file(vcf_file):
    # GIVEN: a good vcf file
    # WHEN-THEN: we validate the file it raises no errors
    validate(vcf_file)


def test_bgzip_good_tbi_file(vcf_file):
    # GIVEN: a good vcf TBI file
    # WHEN-THEN: we validate the file it raises no errors
    validate(vcf_file.with_suffix(".gz.tbi"))


################################################################################
# Error if not a BGZIP file (BgzipMagicNumberError)
################################################################################
def test_bgzip_bad_magic_number():
    # GIVEN: a BPM file
    # WHEN-THEN: we try to validate as a BGZIP file we get an error
    with pytest.raises(BgzipMagicNumberError):
        validate(FakeData._data_path / FakeData._illumina_manifest_file)


################################################################################
# Error if file is truncated (BgzipTruncatedFileError)
################################################################################
def test_truncated_bgzip_file(truncated_vcf_file):
    # GIVEN: a vcf file with a small number of bytes removed from the end
    # WHEN-THEN: we validate the file we get a truncation error
    with pytest.raises(BgzipTruncatedFileError):
        validate(truncated_vcf_file)
