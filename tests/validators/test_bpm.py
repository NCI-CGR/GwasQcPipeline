"""Test BPM validation.

Validation is provided by the Illumina parsing library. Since this is a
binary format, it would be really hard to test individual issues. For now I
am just testing that a good file passes.
"""
import pytest

from cgr_gwas_qc.testing.data import FakeData
from cgr_gwas_qc.validators.bpm import BpmMagicNumberError, validate


################################################################################
# Good BPM file validates
################################################################################
def test_bpm_good_bpm_file():
    # GIVEN: a good bpm file
    # WHEN-THEN: we validate the file it raises no errors
    validate(FakeData._data_path / FakeData._illumina_manifest_file)


################################################################################
# Error if not a BPM file (BPMMagicNumberError)
################################################################################
def test_bpm_bad_magic_number():
    # GIVEN: a good gtc file
    # WHEN-THEN: we validate as a BPM file, then it raises an error
    with pytest.raises(BpmMagicNumberError):
        validate(FakeData._data_path / FakeData._test_gtc)
