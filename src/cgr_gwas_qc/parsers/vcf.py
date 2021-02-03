from contextlib import contextmanager
from typing import Iterable

import pysam


@contextmanager
def open(filename, mode: str = "r"):
    vcf_file = VariantFile(filename, mode)
    try:
        yield vcf_file
    finally:
        vcf_file.close()


class VariantFile(pysam.VariantFile):
    def fetch(
        self,
        contig=None,
        start=None,
        stop=None,
        region=None,
        reopen=False,
        end=None,
        reference=None,
    ) -> Iterable:

        # Ensure correct types if present
        contig = str(contig) if contig else None
        start = int(start) if start else None
        stop = int(stop) if stop else None

        # Sanity checks
        vcf_contigs = self.header.contigs.keys()
        if contig is None or contig in vcf_contigs:
            contig = contig
        elif self._fix_contig(contig) in vcf_contigs:
            contig = self._fix_contig(contig)
        else:
            # Chromosome is not in the VCF
            return []

        if start is not None and start < 0:
            # Negative start position, the smallest allowed value is 0.
            return []

        # Run query
        return super().fetch(contig, start, stop, region, reopen, end, reference)

    @staticmethod
    def _fix_contig(contig: str):
        if contig.startswith("chr"):
            # assume the VCF does not start with "chr"
            return contig.replace("chr", "")

        # assume the VCF does start with "chr"
        return "chr" + contig
