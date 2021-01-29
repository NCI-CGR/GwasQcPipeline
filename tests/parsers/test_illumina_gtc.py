import numpy as np


################################################################################
# There should be 31 entries in the TOC.
################################################################################
def test_length_of_gtc_toc(gtc):
    # GIVEN: A gtc file parsed by the Illumina Genotype calls
    # THEN: The number of TOC entries is the number of records (which I identified as 31)
    assert len(gtc.toc_table) == 31


################################################################################
# Base calls {TOP, FWD, PLUS} should equal the number of loci in the BPM file.
################################################################################
def test_get_top_base_calls(gtc, bpm):
    # GIVEN: A gtc file parsed by the Illumina Genotype calls
    # WHEN: we get the base calls for the TOP strand
    base_calls = list(gtc.get_base_calls())

    # THEN: the number of base_calls is equal to the number of SNPs in the BPM file
    assert len(base_calls) == bpm.num_loci


def test_get_fwd_base_calls(gtc, bpm):
    # GIVEN: A gtc file parsed by the Illumina Genotype calls
    # WHEN: we get the base calls for the FWD strand
    base_calls = list(gtc.get_base_calls_forward_strand(bpm.snps, bpm.source_strands))

    # THEN: the number of base_calls is equal to the number of SNPs in the BPM file
    assert len(base_calls) == bpm.num_loci


def test_get_plus_base_calls(gtc, bpm):
    # GIVEN: A gtc file parsed by the Illumina Genotype calls
    # WHEN: we get the base calls for the PLUS strand
    base_calls = list(gtc.get_base_calls_plus_strand(bpm.snps, bpm.ref_strands))

    # THEN: the number of base_calls is equal to the number of SNPs in the BPM file
    assert len(base_calls) == bpm.num_loci


################################################################################
# Genotype scores should be a numpy array with values between 0 - 1.
################################################################################
def test_genotype_scores(gtc):
    # GIVEN: A gtc file parsed by the Illumina Genotype calls
    # WHEN: we get the genotype scores
    genotype_scores = gtc.get_genotype_scores()

    # THEN:
    # The returned object is a numpy array
    assert isinstance(genotype_scores, np.ndarray)
    # The minimum value is >= 0
    assert 0 <= min(genotype_scores)
    # The minimum value is <= 0
    assert max(genotype_scores) <= 1
