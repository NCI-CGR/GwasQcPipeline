"""Create a list of markers that do not match the VCF."""
from collections import Counter
from pathlib import Path
from textwrap import dedent
from typing import Generator, Optional, Tuple

import pysam
import typer

from cgr_gwas_qc.parsers.illumina import complement

app = typer.Typer(add_completion=False)
counter: Counter = Counter()  # match, remove, flip, bad_chrom, ambiguous, indel, duplicate, missing
unique_snps = set()


@app.command()
def main(
    bim_file: Path = typer.Argument(
        ..., help="A multi-sample plink BIM file.", exists=True, readable=True
    ),
    vcf_file: Path = typer.Argument(..., help="The 1KG VCF file.", exists=True, readable=True),
    snp_removal_list: Path = typer.Argument(
        ..., help="Text file to save the list of makers to remove.", file_okay=True, writable=True,
    ),
    output_bim: Optional[Path] = typer.Argument(
        None,
        help="A BIM file where alleles are corrected to be on the same strand as the VCF.",
        file_okay=True,
        writable=True,
    ),
):
    """Compares a plink BIM file with a VCF and flags SNPs for removal.

    This script creates a filter list that can be used with ``plink --exclude
    <file>``. Filter critera include:

        - SNPs that don't match (and complement doesn't match) the VCF (typically from 1000G)
        - SNPs that are duplicated
        - SNPs that have ambiguous alleles (A/T or C/G)
        - Indels (included on 5M chip)
        - SNPs on chromosomes other than 1-22,X

    This script also outputs a BIM file where flipped alleles are corrected
    to match the VCF.
    """
    if output_bim is None:
        output_bim = bim_file.with_suffix(".vcfStrand.bim")

    vcf = pysam.VariantFile(vcf_file, "r")
    remove = snp_removal_list.open("w")
    bim_out = output_bim.open("w")

    try:
        for row in BimFile(bim_file):
            if no_match(row, vcf):
                remove.write(row.vid_string())
            bim_out.write(row.bim_string())
    finally:
        vcf.close()
        remove.close()
        bim_out.close()

    typer.echo(
        dedent(
            f"""
            ########################################
            #{"Summary":^38}#
            ########################################
            {"Number SNPs Match:": <30}{counter["match"]:>6}
            {"  Number SNPs Flipped:": <30}{counter["flip"]:>6}
            {"Number SNPs to Remove:": <30}{counter["remove"]:>6,}
            {"  Unrecognized Chromosome:": <30}{counter["bad_chrom"]:>6,}
            {"  Ambiguous Alleles:": <30}{counter["ambiguous"]:>6,}
            {"  Insertion or Deletion:": <30}{counter["indel"]:>6,}
            {"  Duplicate SNPs:": <30}{counter["duplicate"]:>6,}
            {"  Not in VCF:": <30}{counter["missing"]:>6,}
            """
        )
    )


def no_match(row, vcf) -> bool:
    """Determine if markers is problematic or not in the VCF."""
    if row.not_major_chrom():
        typer.echo(f"Unrecognized chrom [{row.chrom}]: {row.variant_id}")
        counter["remove"] += 1
        counter["bad_chrom"] += 1
        return True

    if row.is_ambiguous_allele():
        typer.echo(f"Ambiguous alleles [{'/'.join(row.alleles)}]: {row.variant_id}")
        counter["remove"] += 1
        counter["ambiguous"] += 1
        return True

    if row.is_indel():
        typer.echo(f"Indel: {row.variant_id}")
        counter["remove"] += 1
        counter["indel"] += 1
        return True

    for record in vcf.fetch(row.chrom_decode, row.coordinate - 1, row.coordinate):
        if record.id is None:
            record.id = f"{record.chrom}_{record.pos}:{record.ref}:{record.alts[0]}"

        if any("<" in alt for alt in record.alts):
            continue

        if row.coordinate == record.pos and sorted(row.alleles) == sorted(record.alleles):
            row.variant_id = record.id

            if is_duplicate(row):
                counter["remove"] += 1
                counter["duplicate"] += 1
                return True

            counter["match"] += 1
            return False

        if row.coordinate == record.pos and sorted(row.allele_complements) == sorted(
            record.alleles
        ):
            row.variant_id = record.id

            if is_duplicate(row):
                counter["remove"] += 1
                counter["duplicate"] += 1
                return True

            _alleles = row.alleles
            row.set_complements()
            typer.echo(
                f"Flipping alleles [{'/'.join(_alleles)} -> {'/'.join(row.alleles)}]: {row.variant_id}"
            )
            counter["match"] += 1
            counter["flip"] += 1
            return False

    typer.echo(f"Not in VCF: {row.variant_id}")
    counter["remove"] += 1
    counter["missing"] += 1
    return True


def is_duplicate(row: "BimRecord") -> bool:
    if row.variant_id in unique_snps:
        row.variant_id += "_duplicate"
        return True
    unique_snps.add(row.variant_id)
    return False


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
        return complement(self.allele_1), complement(self.allele_2)

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
    app()
