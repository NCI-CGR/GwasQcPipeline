import pytest

from cgr_gwas_qc.parsers import vcf


@pytest.fixture(scope="module")
def vcf_fh(vcf_file):
    return vcf.VcfFile(vcf_file)


@pytest.mark.parametrize("chrom", [1, "1", "chr1", "X", "chrX"])
def test_vcf_file_fetch_different_chrom_formats(vcf_fh, chrom):
    # GIVEN: VCF file handler
    # THEN: I can fetch different chromosome codes and they will be handled
    # gracefully.
    assert vcf_fh.fetch(chrom)


@pytest.mark.parametrize("chrom", [-1, "XY", "NA"])
def test_vcf_file_fetch_nonsense_chrom(vcf_fh, chrom):
    # GIVEN: VCF file handler
    # THEN: I can fetch nonsense chromosome codes and they will be
    # handled gracefully.
    assert [] == list(vcf_fh.fetch(chrom))


@pytest.mark.parametrize(
    "record,result",
    [
        (vcf.VcfRecord(None, "1", 1, "T", ("G", "C")), True),
        (vcf.VcfRecord(None, "1", 1, "T", ("G",)), False),
    ],
)
def test_vcf_record_is_multiallelic(record, result):
    assert record.is_multiallelic() == result


@pytest.mark.parametrize(
    "record,result",
    [
        (vcf.VcfRecord(None, "1", 1, "T", ("G",)), True),
        (vcf.VcfRecord(None, "1", 1, "T", ("GC",)), False),
        (vcf.VcfRecord(None, "1", 1, "TAA", ("C",)), False),
    ],
)
def test_vcf_record_is_snp(record, result):
    assert record.is_snp() == result
