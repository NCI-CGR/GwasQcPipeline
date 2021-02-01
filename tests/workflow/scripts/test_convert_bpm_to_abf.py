"""Test the conversion of BPM to ABF.

For contamination checking we need to have the population level allele
frequencies. We pull these allele frequencies from the 1KG project for each
SNP in the BPM.
"""
import warnings
from math import isclose
from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal
from typer.testing import CliRunner

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

    obs_abf = pd.read_csv(file_out, sep="\t").fillna(0)  # legacy uses 0.0 instead of NA.
    exp_abf = pd.read_csv("tests/data/small_manifest.AF.abf.txt", sep="\t")
    assert_frame_equal(obs_abf, exp_abf)


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
# GitLab Issue #30 failed regression test with real data.
################################################################################
@pytest.mark.real_data
@pytest.mark.regression
def test_bpm2abf_real_data(tmp_path):
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


@pytest.mark.parametrize("pos", [-10, 0], ids=["negative", "zero"])
def test_bpm2abf_impossible_positions(vcf_file, pos):
    from cgr_gwas_qc.parsers.vcf import VariantFile
    from cgr_gwas_qc.workflow.scripts.bpm2abf import Variant, get_abf_from_vcf

    # GIVEN: The 1KG VCF and a variant with and impossible position.
    vcf = VariantFile(vcf_file)
    var = Variant("1", pos, "test_snp", "[A/C]", "+")

    # WHEN: I try to look up the variant in the VCF
    res = get_abf_from_vcf(vcf, "AF", var)

    # THEN: I get no exceptions and I return None
    assert res is None


################################################################################
# GRCh38 Support
################################################################################
@pytest.mark.real_data
def test_bpm2abf_GRCh38_AF(tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data_cache = RealData(GRCh_version=38)

    bpm_file = data_cache / data_cache._illumina_manifest_file
    vcf_file = data_cache / data_cache._thousand_genome_vcf
    file_out = tmp_path / bpm_file.with_suffix(".abf.txt").name
    results = runner.invoke(
        app, [bpm_file.as_posix(), vcf_file.as_posix(), "AF", file_out.as_posix()]
    )
    assert results.exit_code == 0

    obs_abf = pd.read_csv(file_out, sep="\t")
    assert obs_abf.notna().mean().ABF >= 0.90  # Make sure >=90% of SNPs are not missing
