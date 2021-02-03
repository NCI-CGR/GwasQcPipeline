"""Compares a plink BIM file with the 1KG VCF and flags SNPs for removal."""
import re
from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal
from typer.testing import CliRunner

from cgr_gwas_qc.parsers.bim import BimRecord
from cgr_gwas_qc.parsers.vcf import VariantFile
from cgr_gwas_qc.testing import file_hashes_equal

runner = CliRunner()


################################################################################
# Compare to outputs from old script.
################################################################################
@pytest.mark.regression
def test_bim_filter_vcf(bim_file, vcf_file, tmp_path):
    from cgr_gwas_qc.workflow.scripts.bim_filter_vcf import app

    expected_snps_to_remove = Path("tests/data/SnpsToRemove.txt")
    expected_strand_corrected_bim = Path("tests/data/plink/samples.vcfStrand.bim")

    # Run script and create observed outputs
    obs_snps_to_remove = tmp_path / "snps_to_remove.txt"
    obs_strand_corrected_bim = tmp_path / "vcfStrand.bim"

    results = runner.invoke(
        app,
        [
            bim_file.as_posix(),
            vcf_file.as_posix(),
            obs_snps_to_remove.as_posix(),
            obs_strand_corrected_bim.as_posix(),
        ],
    )
    assert results.exit_code == 0

    # Make sure observed and expected outputs are exactly equal
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
    print(results.stdout)
    assert re.match(regexs["match"], results.stdout)
    assert re.match(regexs["flip"], results.stdout)
    assert re.match(regexs["remove"], results.stdout)
    assert re.match(regexs["chrom"], results.stdout)
    assert re.match(regexs["position"], results.stdout)
    assert re.match(regexs["ambiguous"], results.stdout)
    assert re.match(regexs["indel"], results.stdout)
    assert re.match(regexs["duplicate"], results.stdout)
    assert re.match(regexs["missing"], results.stdout)


################################################################################
# Test specific examples
################################################################################
ok_bim_records = [  # record, update_msg
    # Ignore snps not in VCF
    (BimRecord("1", "rs00001", 0, 12341234, "A", "G"), "missing"),
    (BimRecord("1", "rs00001", 0, 12341234, "T", "C,A"), "missing"),
    # Match Snps in VCF
    (BimRecord("1", "rs148369513", 0, 167042622, "C", "T"), "exact_match"),
    (
        BimRecord("1", "rs148369513", 0, 167042622, "T", "C"),
        "exact_match",
    ),  # Match even if the alleles are switched
    (BimRecord("1", "rs148369513", 0, 167042622, "G", "A"), "flip"),  # Match even if the complement
    (
        BimRecord("1", "rs148369513", 0, 167042622, "A", "G"),
        "flip",
    ),  # Match even if the complement and switched
]


@pytest.mark.parametrize("record,update_msg", ok_bim_records)
def test_update_record_with_vcf(vcf_file, record, update_msg, monkeypatch):
    from cgr_gwas_qc.workflow.scripts.bim_filter_vcf import update_bim_record_with_vcf

    monkeypatch.setattr("cgr_gwas_qc.workflow.scripts.bim_filter_vcf.unique_snps", set())

    vcf = VariantFile(vcf_file)
    res = update_bim_record_with_vcf(record, vcf)
    assert update_msg == res


def test_update_record_with_vcf_duplicates(vcf_file, monkeypatch):
    from cgr_gwas_qc.workflow.scripts.bim_filter_vcf import update_bim_record_with_vcf

    monkeypatch.setattr("cgr_gwas_qc.workflow.scripts.bim_filter_vcf.unique_snps", set())

    vcf = VariantFile(vcf_file)
    record1 = BimRecord("1", "rs148369513", 0, 167042622, "C", "T")
    record2 = BimRecord("1", "GSA-rs148369513", 0, 167042622, "C", "T")
    messages = [update_bim_record_with_vcf(record1, vcf), update_bim_record_with_vcf(record2, vcf)]
    assert ["exact_match", "duplicate"] == messages
