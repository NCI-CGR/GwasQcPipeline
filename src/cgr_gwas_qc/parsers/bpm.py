"""Parser for the Illumina BPM format."""
import re
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Generator, List, Optional

from cgr_gwas_qc.parsers import CgrFile
from cgr_gwas_qc.parsers.illumina.IlluminaBeadArrayFiles import (
    BeadPoolManifest,
    RefStrand,
    SourceStrand,
    complement,
)


@contextmanager
def open(filename):
    bpm_data = BpmFile(filename)
    try:
        yield bpm_data
    finally:
        del bpm_data


class BpmFile(CgrFile):
    def __init__(self, filename):
        self.fileobj = BeadPoolManifest(filename)
        self.num_loci = self.fileobj.num_loci
        self.manifest_name = self.fileobj.manifest_name
        self.control_config = self.fileobj.control_config

    def __iter__(self) -> Generator["BpmRecord", None, None]:
        for values in zip(
            self.fileobj.names,
            self.fileobj.snps,
            self.fileobj.chroms,
            self.fileobj.map_infos,
            self.fileobj.ref_strands,
            self.fileobj.source_strands,
        ):
            yield BpmRecord(*values)


@dataclass
class BpmRecord:
    id: str
    snp: str
    chrom: str
    pos: int
    _ref_strand: int
    _source_strand: int
    A_allele: Optional[str] = None
    A_allele_complement: Optional[str] = None
    B_allele: Optional[str] = None
    B_allele_complement: Optional[str] = None
    ref_strand: Optional[str] = None
    source_strand: Optional[str] = None

    def __post_init__(self):
        self.A_allele, self.B_allele = self._split_alleles()
        self.A_allele_complement, self.B_allele_complement = self._get_complements()
        self.ref_strand = RefStrand.to_string(self._ref_strand)
        self.source_strand = SourceStrand.to_string(self._source_strand)

    def is_ambiguous(self) -> bool:
        return sorted([self.A_allele, self.B_allele]) in [["A", "T"], ["C", "G"]]

    def is_indel(self) -> bool:
        return self.A_allele in ["D", "I"] or self.B_allele in ["D", "I"]

    def get_record_problems(self) -> List[str]:
        problems = []
        if self.pos < 1:
            problems.append("bad_position")

        if self.is_ambiguous():
            problems.append("ambiguous_allele")

        if self.is_indel():
            problems.append("indel")

        return problems

    def _split_alleles(self):
        return re.search(r"\[(\w)\/(\w)\]", self.snp).groups()

    def _get_complements(self):
        return complement(self.A_allele), complement(self.B_allele)
