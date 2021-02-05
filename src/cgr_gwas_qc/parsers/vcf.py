from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any, Dict, Iterable, Optional, Tuple

import pysam


@contextmanager
def open(filename, mode: str = "r"):
    vcf_file = VcfFile(filename, mode)
    try:
        yield vcf_file
    finally:
        vcf_file.close()


class VcfFile(pysam.VariantFile):
    def fetch(
        self,
        contig=None,
        start=None,
        stop=None,
        region=None,
        reopen=False,
        end=None,
        reference=None,
    ) -> Iterable["VcfRecord"]:
        """Add better chromosome handling to `pysam.VariantFile.fetch`.

        This overrids `pysam.VariantFile.fetch` to add additional
        sanitization handling before running the query. Also instead of
        returning a `pysam.VariantRecord` I return my own `VcfRecord`
        interface.
        """

        # Ensure correct types if present
        contig = str(contig) if contig else None
        start = int(start) if start else None
        stop = int(stop) if stop else None

        # Try different contig (chromosome) formats. Check the header and see
        # if contig exists and if it doesn't then re-format it and try again.
        # This helps handle ("1" and "chr1").
        vcf_contigs = self.header.contigs.keys()
        if contig is None or contig in vcf_contigs:
            contig = contig
        elif self._fix_contig(contig) in vcf_contigs:
            contig = self._fix_contig(contig)
        else:
            return []  # Chromosome is not in the VCF

        if start is not None and start < 0:
            return []  # Negative start position, the smallest allowed value is 0.

        # Run query
        for record in super().fetch(contig, start, stop, region, reopen, end, reference):
            yield VcfRecord(
                record.id,
                record.chrom,
                record.pos,
                record.ref,
                record.alts,
                record.alleles,
                dict(record.info),
            )

    @staticmethod
    def _fix_contig(contig: str):
        if contig.startswith("chr"):
            # assume the VCF does not start with "chr"
            return contig.replace("chr", "")

        # assume the VCF does start with "chr"
        return "chr" + contig


@dataclass
class VcfRecord:
    id: Optional[str]
    chrom: str
    pos: int
    ref: str
    alts: Tuple[str, ...]
    alleles: Optional[Tuple[str, ...]] = None
    info: Optional[Dict[str, Any]] = None

    def __post_init__(self):
        if self.id is None:
            # If there is no ID then use record information
            _alts = ",".join(self.alts)
            self.id = f"{self.chrom}_{self.pos}:{self.ref}:{_alts}"

        if self.alleles is None:
            self.alleles = (self.ref, *self.alts)

    def is_multiallelic(self):
        return len(self.alts) > 1

    def is_snp(self):
        return len(self.ref) == 1 and len(self.alts[0]) == 1
