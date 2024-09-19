"""Compares a plink BIM file with the 1KG VCF and flags SNPs for removal."""

import re

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from cgr_gwas_qc.parsers import vcf
from cgr_gwas_qc.parsers.bim import BimRecord
from cgr_gwas_qc.testing.comparison import file_hashes_equal
from cgr_gwas_qc.workflow.scripts import bim_filter_vcf


################################################################################
# Compare to outputs from old script.
################################################################################
@pytest.mark.regression
def test_bim_filter_vcf(fake_data_cache, tmp_path, capsys):
    obs_snps_to_remove = tmp_path / "snps_to_remove.txt"
    obs_strand_corrected_bim = tmp_path / "vcfStrand.bim"
    bim_filter_vcf.main(
        fake_data_cache / "plink/samples.bim",
        fake_data_cache / "1KG/small_1KG.vcf.gz",
        obs_snps_to_remove,
        obs_strand_corrected_bim,
    )
    captured = capsys.readouterr()

    # Make sure observed and expected outputs are exactly equal
    expected_snps_to_remove = fake_data_cache / "SnpsToRemove.txt"
    expected_strand_corrected_bim = fake_data_cache / "plink/samples.vcfStrand.bim"

    assert file_hashes_equal(expected_snps_to_remove, obs_snps_to_remove)
    assert_frame_equal(
        pd.read_csv(expected_strand_corrected_bim, delim_whitespace=True),
        pd.read_csv(obs_strand_corrected_bim, delim_whitespace=True),
    )

    # Make sure the log (STDOUT) matches what we expect
    regexs = {
        "match": re.compile(r".*Number SNPs Match:\s+372.*", re.DOTALL),
        "flip": re.compile(r".*Number SNPs Flipped:\s+0.*", re.DOTALL),
        "remove": re.compile(r".*Number SNPs to Remove:\s+230.*", re.DOTALL),
        "chrom": re.compile(r".*Unrecognized Chromosome:\s+8.*", re.DOTALL),
        "position": re.compile(r".*Impossible Position:\s+0.*", re.DOTALL),
        "ambiguous": re.compile(r".*Ambiguous Alleles:\s+38.*", re.DOTALL),
        "indel": re.compile(r".*Insertion or Deletion:\s+156.*", re.DOTALL),
        "duplicate": re.compile(r".*Duplicate SNPs:\s+0.*", re.DOTALL),
        "missing": re.compile(r".*Not in VCF:\s+31.*", re.DOTALL),
    }
    assert re.match(regexs["match"], captured.out)
    assert re.match(regexs["flip"], captured.out)
    assert re.match(regexs["remove"], captured.out)
    assert re.match(regexs["chrom"], captured.out)
    assert re.match(regexs["position"], captured.out)
    assert re.match(regexs["ambiguous"], captured.out)
    assert re.match(regexs["indel"], captured.out)
    assert re.match(regexs["duplicate"], captured.out)
    assert re.match(regexs["missing"], captured.out)


################################################################################
# Test specific examples
################################################################################
@pytest.mark.parametrize(
    "record,update_msg",
    [  # record, update_msg
        # Ignore snps not in VCF
        (BimRecord("rs00001", "1", 12341234, "A", "G"), "missing"),
        (BimRecord("rs00001", "1", 12341234, "T", "C,A"), "missing"),
        # Match Snps in VCF
        (BimRecord("rs148369513", "1", 167042622, "C", "T"), "exact_match"),
        (
            BimRecord("rs148369513", "1", 167042622, "T", "C"),
            "exact_match",
        ),  # Match even if the alleles are switched
        (
            BimRecord("rs148369513", "1", 167042622, "G", "A"),
            "flip",
        ),  # Match even if the complement
        (
            BimRecord("rs148369513", "1", 167042622, "A", "G"),
            "flip",
        ),  # Match even if the complement and switched
    ],
)
def test_update_record_with_vcf(fake_data_cache, record, update_msg, monkeypatch):
    monkeypatch.setattr("cgr_gwas_qc.workflow.scripts.bim_filter_vcf.unique_snps", set())

    with vcf.open(fake_data_cache / "1KG/small_1KG.vcf.gz") as vcf_fh:
        res = bim_filter_vcf.update_bim_record_with_vcf(record, vcf_fh)

    assert update_msg == res


def test_update_record_with_vcf_duplicates(fake_data_cache, monkeypatch):
    monkeypatch.setattr("cgr_gwas_qc.workflow.scripts.bim_filter_vcf.unique_snps", set())

    record1 = BimRecord("rs148369513", "1", 167042622, "C", "T")
    record2 = BimRecord("GSA-rs148369513", "1", 167042622, "C", "T")

    with vcf.open(fake_data_cache / "1KG/small_1KG.vcf.gz") as vcf_fh:
        messages = [
            bim_filter_vcf.update_bim_record_with_vcf(record1, vcf_fh),
            bim_filter_vcf.update_bim_record_with_vcf(record2, vcf_fh),
        ]

    assert ["exact_match", "duplicate"] == messages
