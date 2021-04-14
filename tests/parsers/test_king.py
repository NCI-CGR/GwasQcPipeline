import pytest

from cgr_gwas_qc.parsers import king


@pytest.fixture
def king_kin(tmp_path):
    outfile = tmp_path / "king.kin0"
    outfile.write_text(
        "FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tKinship\n"
        "S0001\tS0001\tS0002\tS0002\t645397\t0.0501\t0.0457\t-0.1305\n"  # unrelated
        "S0002\tS0002\tS0003\tS0003\t645594\t0.0564\t0.0280\t0.36\n"  # identical
        "S0003\tS0003\tS0004\tS0004\t645534\t0.0563\t0.0280\t0.2\n"  # 1st degree
        "S0004\tS0004\tS0005\tS0005\t645534\t0.0563\t0.0280\t0.1\n"  # 2nd degree
        "S0004\tS0004\tS0005\tS0005\t645534\t0.0563\t0.0280\t0.05\n"  # 3rd degree
    )
    return outfile


def test_read_kinship(king_kin):
    df = king.read_kinship(king_kin)
    assert ["UN", "ID", "D1", "D2", "D3"] == df.relationship.tolist()
