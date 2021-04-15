import pytest

from cgr_gwas_qc.parsers import graf
from cgr_gwas_qc.testing.data import FakeData


@pytest.fixture
def graf_relatedness_tsv(tmp_path):
    return FakeData() / "graf/relatedness.tsv"


def test_read_relatedness_sort_ids(graf_relatedness_tsv):
    """Make sure IDs get sorted alphanumerically."""
    df = graf.read_relatedness(graf_relatedness_tsv)
    ids = df.loc[0, ["ID1", "ID2"]].tolist()
    assert ["SP00001", "SP00002"] == ids


def test_read_relatedness_no_spaces_in_columns(graf_relatedness_tsv):
    """Make sure IDs get sorted alphanumerically."""
    df = graf.read_relatedness(graf_relatedness_tsv)
    assert all(" " not in colname for colname in df.columns)


def test_read_relatedness_relationships(graf_relatedness_tsv):
    """Make sure IDs get sorted alphanumerically."""
    df = graf.read_relatedness(graf_relatedness_tsv).set_index(["ID1", "ID2"])

    assert "UN" == df.loc[("SP00001", "SP00002"), "relationship"]
    assert "PO" == df.loc[("SP00001", "SP00003"), "relationship"]
    assert "D2" == df.loc[("SP00001", "SP00004"), "relationship"]
    assert "D3" == df.loc[("SP00001", "SP00005"), "relationship"]
    assert "ID" == df.loc[("SP00002", "SP00003"), "relationship"]
