from typing import Union

import pysam


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
    ):
        # Ensure correct types if present
        contig = str(contig) if contig else None
        start = int(start) if start else None
        stop = int(stop) if stop else None

        try:
            return super().fetch(contig, start, stop, region, reopen, end, reference)
        except ValueError:
            # Mostlikely this is due to a mismatch in chromosome code between
            # `contig` and the VCF. Adjust contig to match VCF.
            contig = self._fix_contig(contig)
            return super().fetch(contig, start, stop, region, reopen, end, reference)

    def contains_contig(self, contig: Union[int, str]):
        contig = str(contig)
        vcf_contigs = self.header.contigs.keys()
        return contig in vcf_contigs or self._fix_contig(contig) in vcf_contigs

    @staticmethod
    def _fix_contig(contig: str):
        if contig.startswith("chr"):
            # assume the VCF does not start with "chr"
            return contig.replace("chr", "")

        # assume the VCF does start with "chr"
        return "chr" + contig
