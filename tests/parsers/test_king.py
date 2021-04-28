from cgr_gwas_qc.parsers import king
from cgr_gwas_qc.testing.data import FakeData


def test_read_related():
    king_file = FakeData() / "king/king.kin0"
    df = king.read_related(king_file).set_index(["ID1", "ID2"])

    assert "PO" == df.loc[("SP00001", "SP00002"), "relationship"]
    assert "D2" == df.loc[("SP00001", "SP00003"), "relationship"]
    assert "ID" == df.loc[("SP00004", "SP00005"), "relationship"]
    assert "D3" == df.loc[("SP00006", "SP00007"), "relationship"]
