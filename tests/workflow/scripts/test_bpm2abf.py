"""Test the conversion of BPM to ABF.

For contamination checking we need to have the population level allele
frequencies. We pull these allele frequencies from the 1KG project for each
SNP in the BPM.
"""

import pandas as pd
import pytest
from typer.testing import CliRunner

from cgr_gwas_qc.parsers import bpm, vcf
from cgr_gwas_qc.workflow.scripts import bpm2abf

runner = CliRunner()


################################################################################
# Error if population is not in the VCF.
################################################################################
@pytest.mark.regression
def test_bpm2abf_unknown_population(fake_data_cache, tmp_path):
    """Check what happens if the given population is not in the VCF

    It will exit and return a list of possible populations.
    """
    test_abf = tmp_path / "small_manifest.abf.txt"
    results = runner.invoke(
        bpm2abf.app,
        [
            (fake_data_cache / "illumina/bpm/small_manifest.bpm").as_posix(),
            (fake_data_cache / "1KG/small_1KG.vcf.gz").as_posix(),
            "ANY_STRING",
            test_abf.as_posix(),
        ],
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
    assert AF == bpm2abf.get_abf_from_vcf(bpm_record, vcf_mock(vcf_records), "AF")


################################################################################
# Compare to outputs from old script.
################################################################################
@pytest.mark.regression
def test_bpm2abf_AF(fake_data_cache, tmp_path):
    test_abf = tmp_path / "small_manifest.abf.txt"
    results = runner.invoke(
        bpm2abf.app,
        [
            (fake_data_cache / "illumina/bpm/small_manifest.bpm").as_posix(),
            (fake_data_cache / "1KG/small_1KG.vcf.gz").as_posix(),
            "AF",
            test_abf.as_posix(),
        ],
    )
    assert results.exit_code == 0

    obs_df = pd.read_csv(test_abf, sep="\t").fillna(0)  # legacy uses 0.0 instead of NA.
    exp_df = pd.read_csv(fake_data_cache / "small_manifest.AF.abf.txt", sep="\t")
    prop_very_different = (abs(obs_df.ABF - exp_df.ABF) > 1e-6).mean()
    assert 0.01 > prop_very_different  # less than 1% are more than 0.000001 different
