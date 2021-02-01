import warnings

import pytest

from cgr_gwas_qc.parsers.vcf import VariantFile
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.real_data
@pytest.fixture(scope="session", params=[37, 38])
def thousand_genomes_vcf(request):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data_cache = RealData(GRCh_version=request.param)

    return VariantFile(data_cache / data_cache._thousand_genome_vcf)


@pytest.mark.real_data
@pytest.mark.parametrize("chrom", [1, "1", "chr1", "X", "chrX"])
def test_vcf_multiple_chrom_codes(thousand_genomes_vcf, chrom):
    # GIVEN: Thousand Genomes VCF
    # THEN: I can fetch with different chromosome codes and they will be
    # handled gracefully.
    assert thousand_genomes_vcf.fetch(chrom, 1_000_000, 1_001_000)


@pytest.mark.real_data
@pytest.mark.parametrize("chrom", [-1, "XY", "NA"])
def test_vcf_nonsense_chrom(thousand_genomes_vcf, chrom):
    # GIVEN: Thousand Genomes VCF
    # THEN: I can fetch with different chromosome codes and they will be
    # handled gracefully.
    with pytest.raises(ValueError):
        thousand_genomes_vcf.fetch(chrom, 1_000_000, 1_001_000)


@pytest.mark.real_data
@pytest.mark.parametrize("chrom", [1, "chr2", "X"])
def test_contains_contig(thousand_genomes_vcf, chrom):
    assert thousand_genomes_vcf.contains_contig(chrom)


@pytest.mark.real_data
@pytest.mark.parametrize("chrom", ["A", 42, "NA"])
def test_not_contains_contig(thousand_genomes_vcf, chrom):
    assert not thousand_genomes_vcf.contains_contig(chrom)
