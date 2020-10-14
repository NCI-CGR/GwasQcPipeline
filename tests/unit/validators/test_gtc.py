"""Test GTC validation.

Validation is provided by the Illumina parsing library. Since this is a
binary format, it would be really hard to test individual issues. For now I
am just testing that a good file passes.
"""
import pytest

from cgr_gwas_qc.validators.gtc import GtcMagicNumberError, GtcTuncatedFileError, validate


def test_good_gtc_file_vaidates(gtc_file):
    validate(gtc_file)


def test_bad_magic_number(bpm_file):
    with pytest.raises(GtcMagicNumberError):
        validate(bpm_file)


@pytest.fixture(params=[1, 2, 4])
def truncated_gtc(tmp_path, gtc_file, request):
    data = gtc_file.read_bytes()
    trunc_file = tmp_path / "truncated.gtc"
    with trunc_file.open("wb") as fh:
        fh.write(data[: -request.param])
    return trunc_file


def test_truncated_file(truncated_gtc):
    with pytest.raises(GtcTuncatedFileError):
        validate(truncated_gtc)
