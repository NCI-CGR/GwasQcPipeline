"""Test the convertion of BPM to ABF.

For contamination checking we need to have the population level allele
frequencies. We pull these allele frequencies from the 1KG project for each
SNP in the BPM.
"""
from pathlib import Path

import pytest
from typer.testing import CliRunner

from cgr_gwas_qc.cli.bpm2abf import app
from cgr_gwas_qc.testing import file_rows_almost_equal

runner = CliRunner()


@pytest.mark.regression
def test_bpm2abf_AF(bpm_file, vcf_file, tmpdir):
    file_out = Path(tmpdir) / bpm_file.with_suffix(".abf.txt").name
    results = runner.invoke(
        app, [bpm_file.as_posix(), vcf_file.as_posix(), "AF", file_out.as_posix()]
    )
    assert results.exit_code == 0

    obs_abf = file_out
    exp_abf = Path("tests/data/small_manifest.AF.abf.txt")
    assert file_rows_almost_equal(obs_abf, exp_abf, fuzzy_col=1, sep="\t", header=True)


@pytest.mark.regression
def test_bpm2abf_unknown_population(bpm_file, vcf_file, tmpdir):
    """Check what happens if the given population is not in the VCF

    It will exit and return a list of possible populations.
    """
    file_out = Path(tmpdir) / bpm_file.with_suffix(".abf.txt").name
    results = runner.invoke(
        app, [bpm_file.as_posix(), vcf_file.as_posix(), "ANY_STRING", file_out.as_posix()]
    )
    assert results.exit_code == 1
    assert "Population must be one of: " in results.stdout
