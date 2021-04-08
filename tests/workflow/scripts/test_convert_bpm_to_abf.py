"""Test the conversion of BPM to ABF.

For contamination checking we need to have the population level allele
frequencies. We pull these allele frequencies from the 1KG project for each
SNP in the BPM.
"""
from pathlib import Path

import pandas as pd
import pytest
from typer.testing import CliRunner

from cgr_gwas_qc.parsers import bpm, vcf
from cgr_gwas_qc.testing.data import FakeData
from cgr_gwas_qc.workflow.scripts.bpm2abf import app

runner = CliRunner()


@pytest.fixture(scope="module")
def vcf_file():
    return FakeData._data_path / FakeData._thousand_genome_vcf


@pytest.fixture(scope="module")
def bpm_file():
    return FakeData._data_path / FakeData._illumina_manifest_file


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
# Check main logic
################################################################################
bpm_vcf_records = [
    pytest.param(  # Perfect Match
        bpm.BpmRecord("rs123", "1", 100, "A", "G"),
        [vcf.VcfRecord("rs124", "1", 100, "A", ("G",), None, {"AF": (0.1,)})],
        0.1,
        id="Perfect Match",
    ),
    pytest.param(
        bpm.BpmRecord("rs123", "1", 100, "A", "G"),
        [vcf.VcfRecord("rs124", "1", 100, "G", ("A",), None, {"AF": (0.1,)})],
        0.9,
        id="Alleles Flipped",
    ),
    pytest.param(
        bpm.BpmRecord("rs123", "1", 100, "A", "G"),
        [vcf.VcfRecord("rs124", "1", 100, "G", ("A", "T"), None, {"AF": (0.1, 0.2)})],
        "NA",
        id="Multi-allelic",
    ),
    pytest.param(
        bpm.BpmRecord("rs123", "1", 100, "A", "G"),
        [vcf.VcfRecord("rs124", "1", 100, "G", ("ACGA",), None, {"AF": (0.1,)})],
        "NA",
        id="Indel",
    ),
    pytest.param(
        bpm.BpmRecord("rs123", "1", 100, "A", "G"),
        [vcf.VcfRecord("rs124", "1", 101, "A", ("G",), None, {"AF": (0.1,)})],
        "NA",
        id="Wrong position",
    ),
    pytest.param(
        bpm.BpmRecord("rs123", "1", 100, "A", "G"),
        [
            vcf.VcfRecord("rs124", "1", 100, "A", ("ACGA",), None, {"AF": (0.25,)}),
            vcf.VcfRecord("rs124", "1", 100, "A", ("G",), None, {"AF": (0.1,)}),
        ],
        0.1,
        id="Indel and Perfect Match",
    ),
    pytest.param(
        bpm.BpmRecord("rs123", "1", 100, "A", "G"),
        [
            vcf.VcfRecord("rs124", "1", 100, "A", ("G",), None, {"AF": (0.1,)}),
            vcf.VcfRecord("rs124", "1", 100, "A", ("G",), None, {"AF": (0.7,)}),
        ],
        0.1,
        id="Two Perfect Matches only use first",
    ),
]


@pytest.mark.parametrize("bpm_record,vcf_records,AF", bpm_vcf_records)
def test_get_abf_from_vcf(bpm_record, vcf_records, AF, vcf_mock):
    from cgr_gwas_qc.workflow.scripts.bpm2abf import get_abf_from_vcf

    assert AF == get_abf_from_vcf(bpm_record, vcf_mock(vcf_records), "AF")


################################################################################
# Compare to outputs from old script.
################################################################################
@pytest.mark.regression
def test_bpm2abf_AF(bpm_file, vcf_file, tmp_path):
    file_out = tmp_path / bpm_file.with_suffix(".abf.txt").name
    results = runner.invoke(
        app, [bpm_file.as_posix(), vcf_file.as_posix(), "AF", file_out.as_posix()]
    )
    assert results.exit_code == 0

    obs_df = pd.read_csv(file_out, sep="\t").fillna(0)  # legacy uses 0.0 instead of NA.
    exp_df = pd.read_csv("tests/data/small_manifest.AF.abf.txt", sep="\t")
    prop_very_different = (abs(obs_df.ABF - exp_df.ABF) > 1e-6).mean()
    assert 0.01 > prop_very_different  # less than 1% are more than 0.000001 different
