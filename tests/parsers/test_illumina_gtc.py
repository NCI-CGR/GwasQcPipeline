import numpy as np
import pytest

from cgr_gwas_qc.parsers.illumina import BeadPoolManifest, GenotypeCalls


@pytest.fixture(scope="module")
def gtc(gtc_file) -> GenotypeCalls:
    return GenotypeCalls(gtc_file)


@pytest.fixture(scope="module")
def bpm(bpm_file) -> BeadPoolManifest:
    return BeadPoolManifest(bpm_file)


################################################################################
# There should be 31 entries in the TOC.
################################################################################
def test_length_of_gtc_toc(gtc):
    assert len(gtc.toc_table) == 31


################################################################################
# Base calls {TOP, FWD, PLUS} should equal the number of loci in the BPM file.
################################################################################
def test_get_top_base_calls(gtc, bpm):
    base_calls = list(gtc.get_base_calls())
    assert len(base_calls) == bpm.num_loci


def test_get_fwd_base_calls(gtc, bpm):
    base_calls = list(gtc.get_base_calls_forward_strand(bpm.snps, bpm.source_strands))
    assert len(base_calls) == bpm.num_loci


def test_get_plus_base_calls(gtc, bpm):
    base_calls = list(gtc.get_base_calls_plus_strand(bpm.snps, bpm.ref_strands))
    assert len(base_calls) == bpm.num_loci


################################################################################
# Genotype scores should be a numpy array with values between 0 - 1.
################################################################################
def test_genotype_scores(gtc):
    genotype_scores = gtc.get_genotype_scores()
    assert isinstance(genotype_scores, np.ndarray)
    assert 0 <= min(genotype_scores)
    assert max(genotype_scores) <= 1
