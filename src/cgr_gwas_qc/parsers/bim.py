"""Parser for the PLINK BIM format."""

from contextlib import contextmanager
from dataclasses import dataclass
from typing import Generator, List, Optional

from .common import CgrBiAllelicVariantRecord, CgrFile


@contextmanager
def open(filename, mode: str = "r"):
    """Note this has to be used as a context manager.

    To open and close while not using a `with` block you must to the
    `BimFile` class directly.
    """
    bim_file = BimFile(filename, mode)
    try:
        yield bim_file
    finally:
        bim_file.close()


class BimFile(CgrFile):
    """Provides an iterable interface to BIM files."""

    endchar: Optional[str] = "\n"
    fields = [
        "encoded_chrom",
        "id",
        "morgans",
        "pos",
        "allele_1",
        "allele_2",
    ]

    def __iter__(self) -> Generator["BimRecord", None, None]:
        for row in self.fileobj:
            data = dict(zip(self.fields, row.strip().split()))

            # Fix types
            data["morgans"] = int(data["morgans"])
            data["pos"] = int(data["pos"])

            data["chrom"] = _decode(data["encoded_chrom"])

            yield BimRecord(**data)


def _decode(chrom):
    """Decodes BIM chromosome code."""
    chrom_codes = {
        "23": "X",
        "24": "Y",
        "25": "XY",
        "26": "MT",
    }
    return chrom_codes.get(chrom, chrom)


@dataclass(eq=False)
class BimRecord(CgrBiAllelicVariantRecord):
    encoded_chrom: Optional[str] = None  # encoded chromosome
    morgans: Optional[int] = None  # position in morgans

    def get_record_problems(self) -> List[str]:
        """Checks the record for common problems.

        A convenience method to check the record for a set of common problems
        and return a list of those problems. Potential problems:
        ["not_major_chrom", "bad_position", "ambiguous_allele", "indel].
        """
        problems = []
        if self.not_major_chrom():
            problems.append("not_major_chrom")

        if self.pos < 1:
            problems.append("bad_position")

        if self.is_ambiguous():
            problems.append("ambiguous_allele")

        if self.is_indel():
            problems.append("indel")

        return problems

    def __str__(self):
        return f"{self.encoded_chrom} {self.id} {self.morgans} {self.pos} {self.allele_1} {self.allele_2}"

    def __eq__(self, other):
        """Compare two BimRecords ignoring the allele order.

        Allele order is not consistent when running PLINK, so I want to
        consider two records equal even if the alleles are swapped.
        """
        return (
            isinstance(self, BimRecord)
            and isinstance(other, BimRecord)
            and self.id == other.id
            and self.chrom == other.chrom
            and self.encoded_chrom == other.encoded_chrom
            and self.morgans == other.morgans
            and self.pos == other.pos
            and (
                (self.allele_1 == other.allele_1 and self.allele_2 == other.allele_2)
                or (self.allele_1 == other.allele_2 and self.allele_2 == other.allele_1)
            )
        )
