"""Test BPM validation.

Validation is provided by the Illumina parsing library. Since this is a
binary format, it would be really hard to test individual issues. For now I
am just testing that a good file passes.
"""
from cgr_gwas_qc.validators.bpm import validate


def test_good_bpm_file_vaidates(bpm_file):
    validate(bpm_file)
