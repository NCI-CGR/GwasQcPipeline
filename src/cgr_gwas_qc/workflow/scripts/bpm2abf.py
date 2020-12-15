#!/usr/bin/env python
"""Pull out the B allele frequency from a given subpopulation in the provided vcf.

B allele frequencies are used when examining sample contaminaton.
"""
import re
from pathlib import Path
from typing import Generator, Optional

import pysam
import typer

from cgr_gwas_qc.parsers.illumina import BeadPoolManifest, complement

app = typer.Typer(add_completion=False)


class Variant:
    def __init__(self, chrom: str, pos: int, name: str, snp: str, strand: int) -> None:
        self.chrom = chrom
        self.pos = pos
        self.name = name
        self.snp = snp
        self.strand = strand
        self.b_allele: Optional[str] = None
        self.parse_snp_to_b_allele()

    def parse_snp_to_b_allele(self):
        """Parse the B allele from the BPM SNP string."""
        try:
            b_allele = re.findall(r"\[\w\/(\w)\]", self.snp)[0]
            self.b_allele = complement(b_allele) if self.strand == 2 else b_allele
        except IndexError:
            typer.echo(f"Problem parsing SNP: {self.snp} @{self.chrom}:{self.pos}")


@app.command()
def main(
    bpm_file: Path = typer.Argument(
        ..., help="Path to the Illumina BPM file.", exists=True, readable=True
    ),
    vcf_file: Path = typer.Argument(
        ..., help="Path to the VCF file to get ABF.", exists=True, readable=True
    ),
    population: str = typer.Argument(..., help="Which population in VCF file to use for ABF."),
    abf_file: Optional[Path] = typer.Argument(
        None, help="Output file name. If none given it will use the `BPM_FILE + .abf.txt`."
    ),
):
    """Pull out the B allele frequency for a given subpopulation.

    Takes variants from the provided `bpm_file` and tries to find their B
    allele frequency in the provided `vcf_file`. The B allele frequencies are
    based on the allele frequency for the provided subpopulation
    `population`. If the provided subpopulation doesn't exist, it will output
    a list of available populations. If the variant does not exists in the
    `vcf_file` then the ABF will be set to 0.0.
    """
    vcf = pysam.VariantFile(vcf_file)
    available_populations = sorted([x for x in vcf.header.info.keys() if "AF" in x])

    if population not in available_populations:
        typer.echo(f"Population must be one of: {', '.join(available_populations)}")
        raise typer.Exit(1)

    if abf_file is None:
        abf_file = bpm_file.with_suffix(".abf.txt")

    with abf_file.open("w") as f:
        f.write("SNP_ID\tABF\n")
        for variant in get_bpm_variants(bpm_file):
            abf = get_abf_from_vcf(vcf, population, variant)
            abf_str = str(abf) if abf is not None else "NA"
            f.write(f"{variant.name}\t{abf_str}\n")
    vcf.close()


def get_bpm_variants(bpm_file: Path) -> Generator[Variant, None, None]:
    """Pull all variants out of the bpm file."""
    bpm = BeadPoolManifest(bpm_file)
    for variant in zip(bpm.chroms, bpm.map_infos, bpm.names, bpm.snps, bpm.source_strands):
        yield Variant(*variant)


def get_abf_from_vcf(vcf: pysam.VariantFile, population: str, variant: Variant) -> Optional[float]:
    """Pull the select popultion B allele frequency from the VCF.

    Tries to find the given variant in the VCF. If the variant exists then it
    returns the B allele frequency for the given population. If the variant
    is multiallelic or dose not exist then it returns `None`.
    """
    if variant.b_allele is None:
        # No B allele
        return None

    if variant.chrom not in vcf.header.contigs.keys():
        # Chromosome not in VCF
        return None

    b_allele = variant.b_allele
    b_allele_c = complement(b_allele)

    for record in vcf.fetch(variant.chrom, variant.pos - 1, variant.pos + 1):
        ref, alts = record.ref, record.alts

        if len(alts) > 1:
            # Skip Multiallelic loci
            continue

        alt = alts[0]

        if variant.pos != record.pos or (
            b_allele not in [ref, alt] and b_allele_c not in [ref, alt]
        ):
            # Different SNP
            continue

        if sorted([ref, alt]) in (["A", "T"], ["C", "G"]):
            # Ambiguous SNP
            continue

        allele_freq = record.info.get(population)[0]

        if b_allele == alt or b_allele_c == alt:
            return allele_freq

        if b_allele == ref or b_allele_c == ref:
            return round(1.0 - allele_freq, 8)

    # No applicable VCF SNP found
    return None


if __name__ == "__main__":
    if "snakemake" in locals():
        main(
            Path(snakemake.input.bpm),  # type: ignore # noqa
            Path(snakemake.input.vcf),  # type: ignore # noqa
            snakemake.params.population,  # type: ignore # noqa
            Path(snakemake.output.abf),  # type: ignore # noqa
        )
    else:
        app()
