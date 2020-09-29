#!/usr/bin/env python
"""Pull out the B allele frequency for a given subpopulation from the provided vcf."""
import re
from functools import partial
from pathlib import Path
from typing import Iterable, Optional, Tuple

import pysam
import typer

from cgr_gwas_qc.parsers.illumina import BeadPoolManifest, complement

app = typer.Typer()


@app.command()
def main(
    bpm_file: Path = typer.Argument(
        ..., help="Path to the Illumina BPM file.", exists=True, readable=True
    ),
    vcf_file: Path = typer.Argument(
        ..., help="Path to the VCF file to get ABF.", exists=True, readable=True
    ),
    population: str = typer.Argument(..., help="Which population in VCF file to use for ABF."),
    abf_file: Optional[Path] = typer.Argument(None, help="Output file name."),
):
    """Pull out the B allele frequency for a given subpopulation.

    Takes variants from the provided `bpm_file` and tries to find their B
    allele frequency in the provided `vcf_file`. The B allele frequencies are
    based on the allele frequency for the provided subpopulation
    `population`. If the variant does not exist in the `VCF`
    """

    vcf = pysam.VariantFile(vcf_file)
    available_populations = sorted([x for x in vcf.header.info.keys() if "AF" in x])

    if population not in available_populations:
        typer.echo(f"Population must be one of: {', '.join(available_populations)}")
        typer.Exit(1)

    if abf_file is None:
        abf_file = bpm_file.with_suffix(".abf.txt")

    with abf_file.open("w") as f:
        f.write("SNP_ID\tABF\n")
        variants = get_bpm_variants(bpm_file)
        for abf_row in map(partial(get_abf_row, vcf, population), variants):
            f.write(abf_row)
    vcf.close()


def get_bpm_variants(bpm_file: Path) -> Iterable[Tuple[str, int, str, str, int]]:
    """Pull all variants out of the bpm file."""
    bpm = BeadPoolManifest(bpm_file)
    return zip(bpm.chroms, bpm.map_infos, bpm.names, bpm.snps, bpm.source_strands)


def get_abf_row(
    vcf: pysam.VariantFile, population: str, chrom: str, pos: int, name: str, snp: str, strand: int
) -> str:
    """Gets the B allele frequency and returns a string with SNP_ID and ABF."""
    if chrom not in vcf.header.contigs.keys():
        return f"{name}\t0.0\n"

    try:
        b_allele = re.findall(r"\[\w\/(\w)\]", snp)[0]
    except IndexError:
        return f"{name}\t0.0\n"

    if strand == 2:
        b_allele = complement(b_allele)

    abf = get_vcf_abf(vcf, population, chrom, pos, b_allele)

    return f"{name}\t{abf}\n"


def get_vcf_abf(
    vcf: pysam.VariantFile, population: str, chrom: str, pos: int, b_allele: str
) -> float:
    """Pull the select popultion B allele frequency from the VCF.

    Tries to find the given variant in the VCF. If the variant exists
    then it returns the B allele frequency for the given population. If the
    variant dose not exist then it returns an ABF of 0.0.
    """
    for record in vcf.fetch(chrom, pos, pos + 1):
        ref, alt = record.ref, record.alts[0]

        if sorted([ref, alt]) in (["A", "T"], ["C", "G"]):
            # ambiguous
            continue

        if b_allele not in [ref, alt] or complement(b_allele) not in [ref, alt]:
            # Different SNP
            continue

        allele_freq = record.info.get(population)[0]

        if b_allele == alt or complement(b_allele) == alt:
            return allele_freq

        if b_allele == ref or complement(b_allele) == ref:
            return 1.0 - allele_freq

    return 0.0


if __name__ == "__main__":
    app()
