#!/usr/bin/env python
"""This script converts from sample GTCs into plink PED files.

MAP File:
    The plink MAP file describes all of the markers (1 per row) and has the
    following fields:

    chrom code: The chromosome name.
    variant_id: The variant identifier.
    position: This is an optional position in (centi)morgans. To ignore we set it to 0.
    coordinate: This is the location in the reference assembly.

PED File:
    The plink PED file describes the genotype calls for a given sample (row) for
    all markers described in the PED file. It has the following fields:

    Family ID: We use the sample ID here.
    Within-family ID: We also use the sample ID here.
    Within-family ID of father: We ignore this field by setting to 0.
    Within-family ID of mother: We ignore this field by setting to 0.
    Sex code: We ignore this field by setting to 0 (unknown).
    Phenotype code: We ignore this field by setting to -9 (missing).
    Allele A marker_1: This is the first allele at marker 1.
    Allele B marker_1: This is the second allele at marker 1.
    ...
    Allele A marker_n: This is the first allele at marker n.
    Allele B marker_n: This is the second allele at marker n.
"""
from pathlib import Path
from typing import Iterator, List, Optional

import typer

from cgr_gwas_qc.parsers.illumina import BeadPoolManifest, GenotypeCalls

app = typer.Typer(add_completion=False)


@app.command()
def main(
    gtc_file: Path = typer.Argument(
        ..., help="Path to a sample's GTC file.", exists=True, readable=True
    ),
    bpm_file: Path = typer.Argument(
        ..., help="Path to the array's BPM manifest.", exists=True, readable=True
    ),
    sample_id: str = typer.Argument(..., help="Sample ID."),
    plink_out_dir: Path = typer.Argument(
        "./", help="Output folder.", metavar="OUTDIR", dir_okay=True
    ),
    strand: str = typer.Option(
        "TOP",
        "--strand",
        "-s",
        help="Which strand to use {TOP, FWD, PLUS}.",
        callback=lambda value: value.lower(),
    ),
    genotype_threshold: Optional[float] = typer.Option(
        None, "--cutoff", help="Genotype call threshold to filter variants (0-1)."
    ),
) -> None:
    """Convert a GTC + BPM into the plink format (PED + MAP)."""
    plink_out_dir.mkdir(exist_ok=True)
    map_file = plink_out_dir / f"{sample_id}.map"
    ped_file = plink_out_dir / f"{sample_id}.ped"

    bpm = BeadPoolManifest(bpm_file)
    gtc = GenotypeCalls(gtc_file)

    save_map(bpm, map_file)
    save_ped(bpm, gtc, sample_id, strand, genotype_threshold, ped_file)


def save_map(bpm: BeadPoolManifest, map_file: Path) -> None:
    """Converts BPM into a plink MAP file."""
    typer.echo(f"Saving MAP file: {map_file}")
    with map_file.open(mode="w") as f:
        for variant in zip(bpm.chroms, bpm.names, bpm.map_infos):
            # See MAP columns above
            f.write("{} {} 0 {}\n".format(*variant))


def save_ped(
    bpm: BeadPoolManifest,
    gtc: GenotypeCalls,
    sample_id: str,
    strand: str,
    threshold: Optional[float],
    ped_file: Path,
) -> None:
    """Converts a BPM and GTC into a plink PED file."""
    typer.echo(f"Saving PED file: {ped_file}")
    with ped_file.open(mode="w") as f:
        if threshold:
            flag_no_genotype = gtc.get_genotype_scores() <= threshold
        else:
            flag_no_genotype = [False] * bpm.num_loci
        base_calls = get_base_calls(bpm, gtc, strand, flag_no_genotype)
        # See PED columns above
        f.write(f"{sample_id} {sample_id} 0 0 0 -9 {' '.join(base_calls)}\n")


def get_base_calls(
    bpm: BeadPoolManifest, gtc: GenotypeCalls, strand: str, flags: List[bool]
) -> Iterator[str]:
    """List of strand specific base calls for each marker.

    Pulls out the per marker base calls using the provided strand.
    """
    if strand == "top":
        base_calls = gtc.get_base_calls()
    elif strand == "fwd":
        base_calls = gtc.get_base_calls_forward_strand(bpm.snps, bpm.source_strands)
    else:
        base_calls = gtc.get_base_calls_plus_strand(bpm.snps, bpm.ref_strands)
    return map(convert_base_call, base_calls, flags)


def convert_base_call(base_call: str, flag_no_call: bool) -> str:
    """Converts base calls to plink format.

    Sets base calls below the genotype score threshold or missing (-) to "0".
    """
    if flag_no_call or base_call in ["--", "-"]:
        return "0 0"

    return "{} {}".format(*base_call)


if __name__ == "__main__":
    app()
