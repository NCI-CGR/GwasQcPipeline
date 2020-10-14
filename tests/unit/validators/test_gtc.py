"""Test GTC validation.

Validation is provided by the Illumina parsing library. Since this is a
binary format, it would be really hard to test individual issues. For now I
am just testing that a good file passes.
"""
import pytest

from cgr_gwas_qc.validators.gtc import GtcMagicNumberError, validate


def test_good_gtc_file_vaidates(gtc_file):
    validate(gtc_file)


def test_bad_magic_number(bpm_file):
    with pytest.raises(GtcMagicNumberError):
        validate(bpm_file)
