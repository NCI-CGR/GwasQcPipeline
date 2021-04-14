import pytest

from cgr_gwas_qc.parsers import graf


@pytest.fixture
def graf_relatedness_tsv(tmp_path):
    outfile = tmp_path / "graf.tsv"

    outfile.write_text(
        "# HG match: number of SNPs with matched genotypes when only homozygous SNPs are counted\n"
        "#\n"
        "# Expected HGMR and AGMR values (%) for different types of relationships:\n"
        "# FS HGMR: 4.9796\n"
        "sample1\tsample2\tsubject1\tsubject2\tsex1\tsex2\tHG match\tHG miss\tHGMR\tAG match\tAG miss\tAGMR\tgeno relation\tped relation\tp_value\n"
        "S0002\tS0001\tS0002\tS0001\t0\t0\t1990\t0\t0.00\t3366\t0\t0.00\tID\tDP\t\n"
        "S0002\tS0003\tS0002\tS0003\t0\t0\t1977\t0\t0.00\t3368\t0\t0.00\tID\tDP\t\n"
        "S0004\tS0005\tS0004\tS0005\t0\t0\t2027\t0\t0.00\t3364\t0\t0.00\tID\tDP\t\n"
    )

    return outfile


def test_read_graf_relatedness_sort_ids(graf_relatedness_tsv):
    """Make sure IDs get sorted alphanumerically."""
    df = graf.read_graf_relatedness(graf_relatedness_tsv)
    ids = df.loc[0, ["ID1", "ID2"]].tolist()
    assert ["S0001", "S0002"] == ids


def test_read_graf_relatedness_no_spaces_in_columns(graf_relatedness_tsv):
    """Make sure IDs get sorted alphanumerically."""
    df = graf.read_graf_relatedness(graf_relatedness_tsv)
    assert all(" " not in colname for colname in df.columns)
