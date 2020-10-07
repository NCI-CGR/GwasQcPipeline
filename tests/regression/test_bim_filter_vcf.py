"""Compares a plink BIM file with the 1KG VCF and flags SNPs for removal."""
import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from cgr_gwas_qc.cli.bim_filter_vcf import app
from cgr_gwas_qc.testing import file_hashes_equal

runner = CliRunner()


@pytest.fixture
def expected_snps_to_remove() -> Path:
    return Path("tests/data/SnpsToRemove.txt")


@pytest.fixture
def expected_strand_corrected_bim() -> Path:
    return Path("tests/data/plink/samples.vcfStrand.bim")


@pytest.fixture
def expected_log() -> Path:
    return Path("tests/data/SnpsToRemove.log")


@pytest.fixture
def simulation_json() -> Path:
    return Path("tests/data/plink/sample.json")


@pytest.mark.regression
def test_bim_filter_vcf(
    bim_file, vcf_file, expected_snps_to_remove, expected_strand_corrected_bim, tmpdir
):
    obs_snps_to_remove = Path(tmpdir) / "snps_to_remove.txt"
    obs_strand_corrected_bim = Path(tmpdir) / "vcfStrand.bim"

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
    assert file_hashes_equal(expected_snps_to_remove, obs_snps_to_remove)
    assert file_hashes_equal(expected_strand_corrected_bim, obs_strand_corrected_bim)

    regexs = {
        "match": re.compile(r".*Number SNPs Match:\s+372.*", re.DOTALL),
        "flip": re.compile(r".*Number SNPs Flipped:\s+0.*", re.DOTALL),
        "remove": re.compile(r".*Number SNPs to Remove:\s+230.*", re.DOTALL),
        "chrom": re.compile(r".*Unrecognized Chromosome:\s+8.*", re.DOTALL),
        "ambiguous": re.compile(r".*Ambiguous Alleles:\s+37.*", re.DOTALL),
        "indel": re.compile(r".*Insertion or Deletion:\s+154.*", re.DOTALL),
        "duplicate": re.compile(r".*Duplicate SNPs:\s+0.*", re.DOTALL),
        "missing": re.compile(r".*Not in VCF:\s+31.*", re.DOTALL),
    }
    assert re.match(regexs["match"], results.stdout)
    assert re.match(regexs["flip"], results.stdout)
    assert re.match(regexs["remove"], results.stdout)
    assert re.match(regexs["chrom"], results.stdout)
    assert re.match(regexs["ambiguous"], results.stdout)
    assert re.match(regexs["indel"], results.stdout)
    assert re.match(regexs["duplicate"], results.stdout)
    assert re.match(regexs["missing"], results.stdout)


bim_rows = [
    ("24 rs00001 0 12341234 A G", True, "bad_chrom"),
    ("0 rs00001 0 12341234 A G", True, "bad_chrom"),
    ("42 rs00001 0 12341234 A G", True, "bad_chrom"),
    ("1 rs00001 0 12341234 A T", True, "ambiguous"),
    ("1 rs00001 0 12341234 T A", True, "ambiguous"),
    ("1 rs00001 0 12341234 G C", True, "ambiguous"),
    ("1 rs00001 0 12341234 C G", True, "ambiguous"),
    ("1 rs00001 0 12341234 I I", True, "indel"),
    ("1 rs00001 0 12341234 0 I", True, "indel"),
    ("1 rs00001 0 12341234 I 0", True, "indel"),
    ("1 rs00001 0 12341234 D D", True, "indel"),
    ("1 rs00001 0 12341234 0 D", True, "indel"),
    ("1 rs00001 0 12341234 D 0", True, "indel"),
    ("1 rs00001 0 12341234 A G", True, "missing"),
    ("1 rs00001 0 12341234 G A", True, "missing"),
    ("1 rs00001 0 12341234 C T", True, "missing"),
    ("1 rs00001 0 12341234 T C", True, "missing"),
    ("1 rs00001 0 12341234 T C,A", True, "missing"),
    ("1 rs148369513 0 167042622 C T", False, "match"),
    ("1 rs148369513 0 167042622 T C", False, "match"),  # Match even if the alleles are switched
    ("1 rs148369513 0 167042622 G A", False, "flip"),  # Match even if the complement
    ("1 rs148369513 0 167042622 A G", False, "flip"),  # Match even if the complement and switched
    ("1 rs148369513 0 167042622 C T", True, "duplicate"),
]


@pytest.mark.parametrize("row,expected_call,status", bim_rows)
def test_no_match(vcf, row, expected_call, status):
    from cgr_gwas_qc.cli.bim_filter_vcf import BimRecord, counter, no_match, unique_snps

    # WARNING: counter is created at import and is shared within the test session.
    # Here I clear it to remove any count artifacts.
    counter.clear()

    bim = BimRecord(row)

    if status == "duplicate":
        # to test duplicate I need to ensure that variant_id is already there.
        unique_snps.add(bim.variant_id)
    else:
        # WARNING: unique_snps is created at import and is shared within the test
        # session. Unless I am testing for duplicates I want to clear it.
        unique_snps.discard(bim.variant_id)

    res = no_match(bim, vcf)

    assert expected_call == res
    assert counter[status] == 1
    assert counter["remove"] == 1 if res else counter["match"] == 1
