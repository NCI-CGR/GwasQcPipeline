import pytest

from cgr_gwas_qc.parsers.illumina import GenotypeCalls


@pytest.fixture(scope="module")
def gtc() -> GenotypeCalls:
    return GenotypeCalls("tests/data/illumina/gtc/small_genotype.gtc")


def test_length_of_gtc_toc(gtc):
    assert len(gtc.toc_table) == 31
