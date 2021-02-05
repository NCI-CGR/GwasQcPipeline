"""Parser for the Illumina BPM format."""
import re
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Generator, List, Optional, Sequence

from cgr_gwas_qc.parsers import CgrBiAllelicVariantRecord, CgrFile
from cgr_gwas_qc.parsers.illumina.IlluminaBeadArrayFiles import (
    BeadPoolManifest,
    RefStrand,
    SourceStrand,
    complement,
)


@contextmanager
def open(filename):
    """Note this has to be used as a context manager.

    To open and close while not using a `with` block you must to the
    `BpmFile` class directly.
    """
    bpm_data = BpmFile(filename)
    try:
        yield bpm_data
    finally:
        del bpm_data


class BpmFile(CgrFile):
    """Provides an iterable interface to BPM files."""

    fields = ["id", "snp", "chrom", "pos", "ref_strand", "source_strand"]

    def __init__(self, filename):
        self.fileobj = BeadPoolManifest(filename)
        self.num_loci = self.fileobj.num_loci
        self.manifest_name = self.fileobj.manifest_name
        self.control_config = self.fileobj.control_config

    def write(self):
        raise NotImplementedError(
            "BPM is a binary format and writing is not currently implemented."
        )

    def __iter__(self) -> Generator["BpmRecord", None, None]:
        for row in zip(
            self.fileobj.names,
            self.fileobj.snps,
            self.fileobj.chroms,
            self.fileobj.map_infos,
            self.fileobj.ref_strands,
            self.fileobj.source_strands,
        ):
            data = dict(zip(self.fields, row))
            data["allele_1"], data["allele_2"] = _split_alleles(data["snp"])
            data["ref_strand"] = RefStrand.to_string(data["ref_strand"])
            data["source_strand"] = SourceStrand.to_string(data["source_strand"])

            yield BpmRecord(**data)


def _split_alleles(snp: str) -> Sequence[str]:
    match = re.search(r"\[(\w)\/(\w)\]", snp)
    return match.groups() if match else ()


@dataclass
class BpmRecord(CgrBiAllelicVariantRecord):
    snp: Optional[str] = None  # in the form [A/B]
    ref_strand: Optional[str] = None  # decoded ref strand [U, +, -], set at init
    source_strand: Optional[str] = None  # decoded source strand [U, F, R], set at init

    @property
    def A_allele(self):
        return self.allele_1

    @property
    def A_allele_complement(self):
        return complement(self.allele_1)

    @property
    def B_allele(self):
        return self.allele_2

    @property
    def B_allele_complement(self):
        return complement(self.allele_2)

    def get_record_problems(self) -> List[str]:
        problems = []
        if self.pos < 1:
            problems.append("bad_position")

        if self.is_ambiguous():
            problems.append("ambiguous_allele")

        if self.is_indel():
            problems.append("indel")

        return problems
