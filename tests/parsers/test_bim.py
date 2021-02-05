from io import StringIO
from textwrap import dedent

import pytest

from cgr_gwas_qc.parsers import bim


def test_bim_reader():
    bim_file = dedent(
        # chrom id morgan coord allele1 allele2
        """\
        1  rs123  0  123  A  G
        1  rs128  0  128  T  G
        """
    )

    expected_records = [
        bim.BimRecord("1", "rs123", 0, 123, "A", "G"),
        bim.BimRecord("1", "rs128", 0, 128, "T", "G"),
    ]

    with bim.open(StringIO(bim_file), "r") as fh:
        records = [record for record in fh]

    assert expected_records == records


def test_bim_writer(tmp_path):
    records = [
        bim.BimRecord("1", "rs123", 0, 123, "A", "G"),
        bim.BimRecord("1", "rs128", 0, 128, "T", "G"),
    ]

    expected = dedent(
        # chrom id morgan coord allele1 allele2
        """\
        1 rs123 0 123 A G
        1 rs128 0 128 T G
        """
    )
    out_file = tmp_path / "test.bim"

    with bim.open(out_file, "w") as fh:
        [fh.write(record) for record in records]

    assert expected == out_file.read_text()


ambiguous = [
    (bim.BimRecord("1", "rs123", 0, 123, "A", "T"), True),
    (bim.BimRecord("1", "rs123", 0, 123, "G", "C"), True),
    (bim.BimRecord("1", "rs123", 0, 123, "A", "C"), False),
    (bim.BimRecord("1", "rs123", 0, 123, "A", "G"), False),
    (bim.BimRecord("1", "rs123", 0, 123, "T", "C"), False),
    (bim.BimRecord("1", "rs123", 0, 123, "T", "G"), False),
]


@pytest.mark.parametrize("record,expected", ambiguous)
def test_is_ambiguous(record, expected):
    assert expected == record.is_ambiguous()


indel = [
    (bim.BimRecord("1", "rs123", 0, 123, "I", "D"), True),
    (bim.BimRecord("1", "rs123", 0, 123, "I", "C"), True),
    (bim.BimRecord("1", "rs123", 0, 123, "A", "D"), True),
    (bim.BimRecord("1", "rs123", 0, 123, "A", "G"), False),
]


@pytest.mark.parametrize("record,expected", indel)
def test_is_indel(record, expected):
    assert expected == record.is_indel()


minor_chrom = [
    (bim.BimRecord("0", "rs123", 0, 123, "A", "G"), True),
    (bim.BimRecord("24", "rs123", 0, 123, "A", "G"), True),
    (bim.BimRecord("25", "rs123", 0, 123, "A", "G"), True),
    (bim.BimRecord("26", "rs123", 0, 123, "A", "G"), True),
    (bim.BimRecord("1", "rs123", 0, 123, "A", "G"), False),
    (bim.BimRecord("23", "rs123", 0, 123, "A", "G"), False),
]


@pytest.mark.parametrize("record,expected", minor_chrom)
def test_not_major_chrom(record, expected):
    assert expected == record.not_major_chrom()


allele_order = [
    (bim.BimRecord("1", "rs123", 0, 123, "A", "T"), ("T", "A")),
    (bim.BimRecord("1", "rs123", 0, 123, "A", "C"), ("C", "A")),
    (bim.BimRecord("1", "rs123", 0, 123, "A", "G"), ("G", "A")),
    (bim.BimRecord("1", "rs123", 0, 123, "T", "C"), ("C", "T")),
    (bim.BimRecord("1", "rs123", 0, 123, "T", "G"), ("G", "T")),
]


@pytest.mark.parametrize("record,expected", allele_order)
def test_switch_allele_order(record, expected):
    record.switch_allele_order()
    assert expected == record.alleles


complements = [
    (bim.BimRecord("1", "rs123", 0, 123, "A", "T"), ("T", "A")),
    (bim.BimRecord("1", "rs123", 0, 123, "A", "C"), ("T", "G")),
    (bim.BimRecord("1", "rs123", 0, 123, "A", "G"), ("T", "C")),
    (bim.BimRecord("1", "rs123", 0, 123, "I", "G"), ("I", "C")),
    (bim.BimRecord("1", "rs123", 0, 123, "D", "G"), ("D", "C")),
]


@pytest.mark.parametrize("record,expected", complements)
def test_complement_alleles(record, expected):
    record.complement_alleles(inplace=True)
    assert expected == record.alleles


problems = [
    (bim.BimRecord("0", "rs123", 0, 123, "A", "C"), ["not_major_chrom"]),
    (bim.BimRecord("1", "rs123", 0, 123, "A", "T"), ["ambiguous_allele"]),
    (bim.BimRecord("1", "rs123", 0, 123, "A", "I"), ["indel"]),
    (bim.BimRecord("1", "rs123", 0, 0, "A", "C"), ["bad_position"]),
    (bim.BimRecord("1", "rs123", 0, 0, "A", "T"), ["ambiguous_allele", "bad_position"]),
]


@pytest.mark.parametrize("record,expected", problems)
def test_get_problems(record, expected):
    problems = record.get_record_problems()
    assert sorted(expected) == sorted(problems)
