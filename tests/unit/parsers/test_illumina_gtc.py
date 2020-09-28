import pytest

from cgr_gwas_qc.parsers.illumina import GenotypeCalls


@pytest.fixture(scope="module")
def gtc(gtc_file) -> GenotypeCalls:
    return GenotypeCalls(gtc_file)


def test_length_of_gtc_toc(gtc):
    assert len(gtc.toc_table) == 31
