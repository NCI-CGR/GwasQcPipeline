"""Parser for the PLINK BIM format."""

from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Generator, List, Optional

from cgr_gwas_qc.parsers import CgrFile
from cgr_gwas_qc.parsers.illumina import complement


@contextmanager
def open(filename, mode: str = "r"):
    bim_file = BimFile(filename, mode)
    try:
        yield bim_file
    finally:
        bim_file.close()


class BimFile(CgrFile):
    endchar: Optional[str] = "\n"

    def __iter__(self) -> Generator["BimRecord", None, None]:
        for row in self.fileobj:
            values = row.strip().split()
            values[2] = int(values[2])  # convert morgans to int
            values[3] = int(values[3])  # convert pos to int
            yield BimRecord(*values)


@dataclass
class BimRecord:
    _chrom: str  # encoded chromosome
    id: str
    morgans: int  # position in morgans
    pos: int  # position in base-coordinate
    allele_1: str
    allele_2: str
    chrom: Optional[str] = None
    alleles: List[str] = field(default_factory=list)

    def __post_init__(self):
        self.chrom = self._decode(self._chrom)
        self.alleles = self._alleles()

    def is_ambiguous(self) -> bool:
        return sorted(self.alleles) in [["A", "T"], ["C", "G"]]

    def is_indel(self) -> bool:
        return any(allele in ["D", "I"] for allele in self.alleles)

    def not_major_chrom(self) -> bool:
        return self._chrom not in [str(x) for x in range(1, 24)]

    def switch_allele_order(self):
        self.allele_1, self.allele_2 = self.allele_2, self.allele_1
        self.alleles = self._alleles()

    def complement_alleles(self, inplace=False):
        _a1 = self._complement(self.allele_1)
        _a2 = self._complement(self.allele_2)
        if inplace:
            self.allele_1 = _a1
            self.allele_2 = _a2
            self.alleles = self._alleles()
        else:
            return _a1, _a2

    def get_record_problems(self) -> List[str]:
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

    def _alleles(self):
        return self.allele_1, self.allele_2

    @staticmethod
    def _decode(chrom):
        chrom_codes = {
            "23": "X",
            "24": "Y",
            "25": "XY",
            "26": "MT",
        }
        return chrom_codes.get(chrom, chrom)

    @staticmethod
    def _complement(allele):
        try:
            return complement(allele)
        except ValueError:
            return allele

    def __str__(self):
        return f"{self._chrom} {self.id} {self.morgans} {self.pos} {self.allele_1} {self.allele_2}"
