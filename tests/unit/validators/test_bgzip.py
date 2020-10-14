"""Test BGZIP validation.

Since this is a binary format, it would be really hard to test individual
issues. For now I am just testing that a good file passes.
"""
import pytest

from cgr_gwas_qc.validators.bgzip import BgzipMagicNumberError, BgzipTruncatedFileError, validate


################################################################################
# Good VCF and TBI files validate ok
################################################################################
def test_good_vcf_file_vaidates(vcf_file):
    validate(vcf_file)


def test_good_tbi_file_vaidates(vcf_file):
    validate(vcf_file.with_suffix(".gz.tbi"))


################################################################################
# Error if not a BGZIP file (BgzipMagicNumberError)
################################################################################
def test_bad_magic_number(bpm_file):
    with pytest.raises(BgzipMagicNumberError):
        validate(bpm_file)


################################################################################
# Error if file is truncated (BgzipTruncatedFileError)
################################################################################
@pytest.fixture(params=[1, 2, 4])
def truncated_vcf(tmp_path, vcf_file, request):
    data = vcf_file.read_bytes()
    trunc_file = tmp_path / "truncated.vcf.gz"
    with trunc_file.open("wb") as fh:
        fh.write(data[: -request.param])
    return trunc_file


def test_truncated_file(truncated_vcf):
    with pytest.raises(BgzipTruncatedFileError):
        validate(truncated_vcf)
