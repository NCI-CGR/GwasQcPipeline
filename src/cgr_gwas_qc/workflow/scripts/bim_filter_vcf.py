#!/usr/bin/env python
"""Create a list of markers that do not match the VCF."""
from collections import Counter
from pathlib import Path
from textwrap import dedent
from typing import List, Optional

import typer

from cgr_gwas_qc.parsers import bim, vcf

app = typer.Typer(add_completion=False)
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
    <file>``. Filter criteria include:

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

    # match, remove, flip, bad_chrom, bad_position, ambiguous, indel, duplicate, missing
    counter: Counter = Counter()

    with vcf.open(vcf_file, "r") as vcf_fh, bim.open(bim_file) as bim_in, bim.open(
        output_bim, "w"
    ) as bim_out, snp_removal_list.open("w") as to_remove:
        for record in bim_in:
            problems = record.get_record_problems()
            if problems:
                counter.update(Counter(problems))
                counter["remove"] += 1
                to_remove.write(record.id + "\n")
            else:
                res = update_bim_record_with_vcf(record, vcf_fh)
                counter[res] += 1

                if res in ["missing", "duplicate"]:
                    counter["remove"] += 1
                    to_remove.write(record.id + "\n")

            bim_out.write(record)

    counter["match"] = counter["exact_match"] + counter["flip"]
    typer.echo(
        dedent(
            f"""
            ########################################
            #{"Summary":^38}#
            ########################################
            {"Number SNPs Match:": <30}{counter["match"]:>6}
            {"  Number SNPs Flipped:": <30}{counter["flip"]:>6}
            {"Number SNPs to Remove:": <30}{counter["remove"]:>6,}
            {"  Unrecognized Chromosome:": <30}{counter["not_major_chrom"]:>6,}
            {"  Impossible Position:": <30}{counter["bad_position"]:>6,}
            {"  Ambiguous Alleles:": <30}{counter["ambiguous_allele"]:>6,}
            {"  Insertion or Deletion:": <30}{counter["indel"]:>6,}
            {"  Duplicate SNPs:": <30}{counter["duplicate"]:>6,}
            {"  Not in VCF:": <30}{counter["missing"]:>6,}
            """
        )
    )


def update_bim_record_with_vcf(b_record: bim.BimRecord, vcf_fh: vcf.VariantFile) -> str:
    for v_record in vcf_fh.fetch(b_record.chrom, b_record.pos - 1, b_record.pos):
        if b_record.pos != v_record.pos:
            # positions aren't the same, this should never happen b/c we are using fetch
            continue

        if len(v_record.alts) > 1:
            # Skip Multiallelic loci
            continue

        if len(v_record.ref) > 1 or len(v_record.alts[0]) > 1:
            # only consider SNVs
            continue

        if v_record.id is None:
            v_record.id = f"{v_record.chrom}_{v_record.pos}:{v_record.ref}:{v_record.alts[0]}"

        # Perfect match
        if alleles_equal(b_record.alleles, v_record.alleles):
            b_record.id = v_record.id
            return "exact_match" if not is_duplicate(b_record) else "duplicate"

        # Complements of alleles match the VCF
        if alleles_equal(b_record.complement_alleles(), v_record.alleles):
            b_record.id = v_record.id
            b_record.complement_alleles(inplace=True)
            return "flip" if not is_duplicate(b_record) else "duplicate"

    return "missing"


def alleles_equal(bim: List[str], vcf: List[str]) -> bool:
    return sorted(bim) == sorted(vcf)


def is_duplicate(record: bim.BimRecord) -> bool:
    if record.id in unique_snps:
        record.id += "_duplicate"
        return True
    unique_snps.add(record.id)
    return False


if __name__ == "__main__":
    if "snakemake" in locals():
        from contextlib import redirect_stdout

        with open(str(snakemake.log), "w") as log:  # type: ignore # noqa
            with redirect_stdout(log):
                main(
                    Path(snakemake.input.bim),  # type: ignore # noqa
                    Path(snakemake.input.vcf),  # type: ignore # noqa
                    Path(snakemake.output.snps_to_remove),  # type: ignore # noqa
                    Path(snakemake.output.bim),  # type: ignore # noqa
                )
    else:
        app()
