"""Test the conversion of BPM to ABF.

For contamination checking we need to have the population level allele
frequencies. We pull these allele frequencies from the 1KG project for each
SNP in the BPM.
"""
from math import isclose
from pathlib import Path

import pytest
from typer.testing import CliRunner

from cgr_gwas_qc.testing import file_rows_almost_equal
from cgr_gwas_qc.testing.data import RealData
from cgr_gwas_qc.workflow.scripts.bpm2abf import app

runner = CliRunner()


################################################################################
# Compare to outputs from old script.
################################################################################
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


################################################################################
# Error if population is not in the VCF.
################################################################################
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


################################################################################
# Issue #30 failed regression test with real data.
################################################################################
@pytest.mark.real_data
@pytest.mark.regression
def test_issue_30(tmp_path):
    # GIVEN: A real Illumina BPM and 1KG vcf file
    data_store = RealData()
    bpm_file = (data_store / data_store._illumina_manifest_file).as_posix()
    vcf_file = (data_store / data_store._thousand_genome_vcf).as_posix()
    obs_abf = tmp_path / "test.abf.txt"

    # WHEN: I run the script using the allele frequencies from "AF"
    results = runner.invoke(
        app, [bpm_file, vcf_file, "AF", obs_abf.as_posix()], catch_exceptions=False
    )
    assert results.exit_code == 0

    # THEN: I expect the results to match the legacy workflow closely
    exp_abf = data_store / "production_outputs/GSAMD-24v1-0_20011747_A1.AF.abf.txt"
    with obs_abf.open() as obs, exp_abf.open() as exp_:
        for i, (obs_row, exp_row) in enumerate(zip(obs, exp_)):
            if i == 0:
                continue  # skip header

            obs_name, obs_value = obs_row.strip().split("\t")
            exp_name, exp_value = exp_row.strip().split("\t")

            assert obs_name == exp_name

            # except when the legacy value is 0.0 due to bugs outlined in Issue #30
            if float(exp_value) == 0.0:
                continue

            assert isclose(float(obs_value), float(exp_value), rel_tol=1e-4)
