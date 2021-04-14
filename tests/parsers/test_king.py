from cgr_gwas_qc.parsers import king
from cgr_gwas_qc.testing.data import FakeData


def test_read_kinship():
    king_file = FakeData() / "king/king.kin0"
    df = king.read_kinship(king_file).set_index(["ID1", "ID2"])

    assert "UN" == df.loc[("SP00001", "SP00002"), "relationship"]
    assert "D1" == df.loc[("SP00001", "SP00003"), "relationship"]
    assert "D2" == df.loc[("SP00001", "SP00004"), "relationship"]
    assert "D3" == df.loc[("SP00001", "SP00005"), "relationship"]
    assert "ID" == df.loc[("SP00002", "SP00003"), "relationship"]
