#!/usr/bin/env python
"""Create a list of markers that do not match the VCF."""
import re
import shutil
from pathlib import Path
from typing import Generator, Optional, Tuple

import typer

from cgr_gwas_qc.parsers.illumina import complement
from cgr_gwas_qc.parsers.vcf import VariantFile

app = typer.Typer(add_completion=False)


@app.command()
def main(
    bed_in: Path = typer.Argument(
        ..., help="A multi-sample plink bed file.", exists=True, readable=True
    ),
    bim_in: Path = typer.Argument(
        ..., help="A multi-sample plink BIM file.", exists=True, readable=True
    ),
    fam_in: Path = typer.Argument(
        ..., help="A multi-sample plink FAM file.", exists=True, readable=True
    ),
    vcf_in: Path = typer.Argument(..., help="The 1KG VCF file.", exists=True, readable=True),
    bed_out: Path = typer.Argument(
        ...,
        help="Just a copy of `bed_in`, done to keep naming schema simple.",
        file_okay=True,
        writable=True,
    ),
    bim_out: Path = typer.Argument(
        ...,
        help="A BIM file where alleles are corrected to be on the same strand as the VCF.",
        file_okay=True,
        writable=True,
    ),
    fam_out: Path = typer.Argument(
        None,
        help="Just a copy of `fam_in`, done to keep naming schema simple.",
        file_okay=True,
        writable=True,
    ),
):
    """Compares a plink BIM file with a VCF and updates IDs to match VCF.

    This script will also use the nucleotide complements of alleles to match
    the VCF.

    .. note::
        The BED_OUT and FAM_OUT are just copies of the input files. PLINK and
        GRAF operate on the entire set of plink binary files so this will
        make using this output easier.

    """
    with VariantFile(vcf_in, "r") as vcf, bim_out.open("w") as bout:
        for row in BimFile(bim_in):
            normalize_variant_id(row)
            update_variant_id_with_vcf(row, vcf)
            bout.write(row.bim_string())

    shutil.copyfile(bed_in, bed_out)
    shutil.copyfile(fam_in, fam_out)


def normalize_variant_id(row: "BimRecord"):
    """Normalize variant id to rsID.

    Sometimes variant IDs are rsID or have an embedded rsID (i.e.,
    GSA-rs1234). If there is some form of rsID in the name go ahead and use
    it as the variants ID.
    """
    rsID = extract_rsID(row.variant_id)
    if rsID:
        row.variant_id = rsID


def extract_rsID(variant: str) -> Optional[str]:
    match = re.search(r"rs\d+", variant)
    return match.group() if match else None


def update_variant_id_with_vcf(row: "BimRecord", vcf: VariantFile) -> None:
    """Update the variant ID using the VCF IDs if present."""
    if row.not_major_chrom() or not vcf.contains_contig(row.chrom_decode):
        # Unrecognized chromosome
        return None

    if row.is_ambiguous_allele():
        # Ambiguous allele
        return None

    if row.is_indel():
        # Indel
        return None

    for record in vcf.fetch(row.chrom_decode, row.coordinate - 1, row.coordinate):
        if record.id is None or extract_rsID(record.id) is None:
            # No rsID to update with
            continue

        if any("<" in alt for alt in record.alts):
            continue

        # Perfect match
        if row.coordinate == record.pos and sorted(row.alleles) == sorted(record.alleles):
            row.variant_id = record.id
            return None

        # Complements of alleles match the VCF
        if row.coordinate == record.pos and sorted(row.allele_complements) == sorted(
            record.alleles
        ):
            row.variant_id = record.id
            row.set_complements()  # Fix the alleles by taking the complement to match VCF record
            return None

    # Variant is not in the VCF
    return None


class BimRecord:
    def __init__(self, row):
        chrom, variant_id, position, coordinate, allele_1, allele_2 = row.strip().split()
        self.chrom: int = int(chrom)
        self.variant_id: str = variant_id
        self.position: int = int(position)
        self.coordinate: int = int(coordinate)
        self.allele_1: str = allele_1
        self.allele_2: str = allele_2
        self.msg: Optional[str] = None

    @property
    def alleles(self) -> Tuple[str, str]:
        return self.allele_1, self.allele_2

    @property
    def allele_complements(self) -> Tuple[str, str]:
        def _complement(allele):
            """Illumina's complement does not work when allele is "0"."""
            try:
                return complement(allele)
            except ValueError:
                return allele

        return _complement(self.allele_1), _complement(self.allele_2)

    @property
    def chrom_decode(self):
        return str(self.chrom) if self.chrom < 23 else "X"

    def bim_string(self) -> str:
        """Convert to a BIM row string ready to write."""
        return "{}\t{}\t{}\t{}\t{}\t{}\n".format(
            self.chrom,
            self.variant_id,
            self.position,
            self.coordinate,
            self.allele_1,
            self.allele_2,
        )

    def vid_string(self) -> str:
        """Create a variant ID string ready to write."""
        return f"{self.variant_id}\n"

    def is_x_chrom(self) -> bool:
        return self.chrom == 23

    def not_major_chrom(self) -> bool:
        return self.chrom == 0 or self.chrom >= 24

    def is_ambiguous_allele(self) -> bool:
        return sorted(self.alleles) in [["A", "T"], ["C", "G"]]

    def is_indel(self) -> bool:
        return any(allele in ["D", "I"] for allele in self.alleles)

    def set_complements(self) -> None:
        """Sets allele_1 and allele_2 to their compliments."""
        self.allele_1, self.allele_2 = self.allele_complements


class BimFile:
    def __init__(self, file_name: Path):
        self.file_name = file_name
        self.fh = file_name.open("r")

    def __iter__(self) -> Generator[BimRecord, None, None]:
        for row in self.fh:
            yield BimRecord(row)


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({k + "_in": Path(v) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({k + "_out": Path(v) for k, v in snakemake.output.items()})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
