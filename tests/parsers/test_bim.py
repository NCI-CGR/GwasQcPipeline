from io import StringIO
from textwrap import dedent

import pytest

from cgr_gwas_qc.parsers import bim


@pytest.fixture(scope="module")
def bim_file():
    return dedent(
        # chrom id morgan coord allele1 allele2
        """\
        1 rs123 0 123 A G
        1 rs128 0 128 T G
        23 rs129 0 129 T G
        """
    )


@pytest.fixture(scope="module")
def bim_records():
    return [
        bim.BimRecord("rs123", "1", 123, "A", "G", "1", 0),
        bim.BimRecord("rs128", "1", 128, "T", "G", "1", 0),
        bim.BimRecord("rs129", "X", 129, "T", "G", "23", 0),
    ]


def test_bim_reader(bim_file, bim_records):
    with bim.open(StringIO(bim_file), "r") as fh:
        records = [record for record in fh]

    assert bim_records == records


def test_bim_writer(bim_file, bim_records, tmp_path):
    out_file = tmp_path / "test.bim"

    with bim.open(out_file, "w") as fh:
        [fh.write(record) for record in bim_records]

    assert bim_file == out_file.read_text()


ambiguous = [
    (bim.BimRecord("rs123", "1", 123, "A", "T"), True),
    (bim.BimRecord("rs123", "1", 123, "G", "C"), True),
    (bim.BimRecord("rs123", "1", 123, "A", "C"), False),
    (bim.BimRecord("rs123", "1", 123, "A", "G"), False),
    (bim.BimRecord("rs123", "1", 123, "T", "C"), False),
    (bim.BimRecord("rs123", "1", 123, "T", "G"), False),
]


@pytest.mark.parametrize("record,expected", ambiguous)
def test_is_ambiguous(record, expected):
    assert expected == record.is_ambiguous()


indel = [
    (bim.BimRecord("rs123", "1", 123, "I", "D"), True),
    (bim.BimRecord("rs123", "1", 123, "I", "C"), True),
    (bim.BimRecord("rs123", "1", 123, "A", "D"), True),
    (bim.BimRecord("rs123", "1", 123, "A", "G"), False),
]


@pytest.mark.parametrize("record,expected", indel)
def test_is_indel(record, expected):
    assert expected == record.is_indel()


minor_chrom = [
    (bim.BimRecord("rs123", "0", 123, "A", "G"), True),
    (bim.BimRecord("rs123", "24", 123, "A", "G"), True),
    (bim.BimRecord("rs123", "25", 123, "A", "G"), True),
    (bim.BimRecord("rs123", "26", 123, "A", "G"), True),
    (bim.BimRecord("rs123", "1", 123, "A", "G"), False),
    (bim.BimRecord("rs123", "X", 123, "A", "G"), False),
]


@pytest.mark.parametrize("record,expected", minor_chrom)
def test_not_major_chrom(record, expected):
    assert expected == record.not_major_chrom()


allele_order = [
    (bim.BimRecord("rs123", "1", 123, "A", "T"), ("T", "A")),
    (bim.BimRecord("rs123", "1", 123, "A", "C"), ("C", "A")),
    (bim.BimRecord("rs123", "1", 123, "A", "G"), ("G", "A")),
    (bim.BimRecord("rs123", "1", 123, "T", "C"), ("C", "T")),
    (bim.BimRecord("rs123", "1", 123, "T", "G"), ("G", "T")),
]


@pytest.mark.parametrize("record,expected", allele_order)
def test_switch_allele_order(record, expected):
    record.switch_allele_order()
    assert expected == record.alleles


complements = [
    (bim.BimRecord("rs123", "1", 123, "A", "T"), ("T", "A")),
    (bim.BimRecord("rs123", "1", 123, "A", "C"), ("T", "G")),
    (bim.BimRecord("rs123", "1", 123, "A", "G"), ("T", "C")),
    (bim.BimRecord("rs123", "1", 123, "I", "G"), ("I", "C")),
    (bim.BimRecord("rs123", "1", 123, "D", "G"), ("D", "C")),
]


@pytest.mark.parametrize("record,expected", complements)
def test_complement_alleles(record, expected):
    record.complement_alleles(inplace=True)
    assert expected == record.alleles


problems = [
    (bim.BimRecord("rs123", "0", 123, "A", "C"), ["not_major_chrom"]),
    (bim.BimRecord("rs123", "1", 123, "A", "T"), ["ambiguous_allele"]),
    (bim.BimRecord("rs123", "1", 123, "A", "I"), ["indel"]),
    (bim.BimRecord("rs123", "1", 0, "A", "C"), ["bad_position"]),
    (bim.BimRecord("rs123", "1", 0, "A", "T"), ["ambiguous_allele", "bad_position"]),
]


@pytest.mark.parametrize("record,expected", problems)
def test_get_problems(record, expected):
    problems = record.get_record_problems()
    assert sorted(expected) == sorted(problems)
