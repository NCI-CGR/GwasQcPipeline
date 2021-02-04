#!/usr/bin/env python
"""Pull out the B allele frequency from a given subpopulation in the provided vcf.

B allele frequencies are used when examining sample contaminaton.
"""
from pathlib import Path
from typing import Optional, Union

import typer

from cgr_gwas_qc.parsers import bpm, vcf

app = typer.Typer(add_completion=False)


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
    `vcf_file` then the ABF will be set to "NA".
    """

    if abf_file is None:
        abf_file = bpm_file.with_suffix(".abf.txt")

    with vcf.open(vcf_file) as vcf_fh:
        available_populations = sorted([x for x in vcf_fh.header.info.keys() if "AF" in x])
        if population not in available_populations:
            typer.echo(f"Population must be one of: {', '.join(available_populations)}")
            raise typer.Exit(1)

        with bpm.open(bpm_file) as bpm_fh, abf_file.open("w") as abf_fh:
            abf_fh.write("SNP_ID\tABF\n")
            for record in bpm_fh:
                b_allele_freq = get_abf_from_vcf(record, vcf_fh, population)
                abf_fh.write(f"{record.id}\t{b_allele_freq}\n")


def get_abf_from_vcf(
    b_record: bpm.BpmRecord, vcf_fh: vcf.VariantFile, population: str
) -> Union[float, str]:
    """Pull the select popultion B allele frequency from the VCF.

    Tries to find the given variant in the VCF. If the variant exists then it
    returns the B allele frequency for the given population. If the variant
    is multiallelic or dose not exist then it returns "NA".
    """
    if b_record.get_record_problems():
        # Problematic record: invalidate coordinate, ambiguous allele, indels
        return "NA"

    for v_record in vcf_fh.fetch(b_record.chrom, b_record.pos - 1, b_record.pos):
        if v_record.pos != b_record.pos:
            # positions aren't the same, this should never happen b/c we are using fetch
            continue

        if len(v_record.alts) > 1:
            # Skip Multiallelic loci
            continue

        ref, alt = v_record.ref, v_record.alts[0]
        if len(ref) > 1 or len(alt) > 1:
            # only consider SNVs
            continue

        allele_freq = v_record.info.get(population)[0]

        if b_record.B_allele == alt or b_record.B_allele_complement == alt:
            return allele_freq

        if b_record.B_allele == ref or b_record.B_allele_complement == ref:
            return round(1.0 - allele_freq, 8)

    # No applicable VCF SNP found
    return "NA"


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({k: Path(v) for k, v in snakemake.output.items()})  # type: ignore # noqa
        defaults.update({k: v for k, v in snakemake.params.items()})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
