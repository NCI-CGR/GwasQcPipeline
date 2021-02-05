import builtins
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional, Tuple

from cgr_gwas_qc.parsers.illumina import complement


class CgrFile(ABC):
    """Base class for interacting with files.

    This is an abstract class for interacting with common genomic file
    formats. It implements open/close/write. And stubs out an iterator to be
    defined in subclasses.
    """

    endchar: Optional[str] = None

    def __init__(self, filename, mode="r"):

        if hasattr(filename, "readlines"):
            # Need for handling StringIO or file objects instead of just file paths.
            self.fileobj = filename
        else:
            self.fileobj = builtins.open(filename, mode)

    def close(self):
        self.fileobj.close()

    def write(self, record, endchar=None):
        """This version of write adds `self.endchar` when writing.

        If the subclass sets `self.endchar` this method will add that
        character when writing record. Useful for automatically adding a
        newline character. This can be directly overridden by setting
        `enchar` to None upon calling.
        """
        payload = str(record)
        _endchar = endchar or self.endchar

        if _endchar is not None:
            payload += _endchar

        self.fileobj.write(payload)

    @abstractmethod
    def __iter__(self):
        raise NotImplementedError


@dataclass
class CgrBiAllelicVariantRecord:
    id: str
    chrom: str
    pos: int
    allele_1: str
    allele_2: str

    @property
    def alleles(self):
        return self.allele_1, self.allele_2

    def is_ambiguous(self) -> bool:
        """Tests if the record's alleles are ambiguous."""
        return sorted(self.alleles) in [["A", "T"], ["C", "G"]]

    def is_indel(self) -> bool:
        """Tests if either of the record's alleles are indels."""
        return any(allele in ["D", "I"] for allele in self.alleles)

    def not_major_chrom(self) -> bool:
        """Tests if the record's chromosome is not 1-22 or X."""
        return self.chrom not in [str(x) for x in range(1, 23)] and self.chrom != "X"

    def switch_allele_order(self) -> None:
        """Updates record switching allele_1 and allele_2"""
        self.allele_1, self.allele_2 = self.allele_2, self.allele_1
        return None

    def complement_alleles(self, inplace=False) -> Optional[Tuple[str, str]]:
        """Gets the complement of allele_1 and allele_2.

        Args:
            inplace: If True then update `allele_1` and `allele_2` with their
              complements. If False then return the complements of `allele_1`
              and `allele_2` but do not update the record. [default False]
        """
        _a1 = self._complement(self.allele_1)
        _a2 = self._complement(self.allele_2)
        if inplace:
            self.allele_1 = _a1
            self.allele_2 = _a2
            return None

        return _a1, _a2

    @staticmethod
    def _complement(allele):
        try:
            return complement(allele)
        except ValueError:
            return allele

    def get_record_problems(self):
        raise NotImplementedError
