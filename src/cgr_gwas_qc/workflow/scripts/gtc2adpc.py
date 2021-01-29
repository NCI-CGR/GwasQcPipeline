"""Converts a GCT + BIM file into Illumina's adpc.bin format."""
from pathlib import Path
from typing import Generator, List, Optional

import typer
from more_itertools import unzip

from cgr_gwas_qc.parsers.illumina import AdpcRecord, AdpcWriter, BeadPoolManifest, GenotypeCalls

app = typer.Typer(add_completion=False)


@app.command()
def main(
    gtc_file: Path = typer.Argument(
        ..., help="Path to a sample(s) GTC file.", exists=True, readable=True
    ),
    bpm_file: Path = typer.Argument(
        ..., help="Path to the array's BPM manifest.", exists=True, readable=True
    ),
    adpc_file: Path = typer.Argument(
        ..., help="Path to output the adpc.bin file.", file_okay=True, writable=True
    ),
    snp_count_file: Optional[Path] = typer.Argument(
        None,
        help="Path to save the number of snps in the BPM manifest.",
        file_okay=True,
        writable=True,
    ),
) -> None:
    """Converts a GCT + BIM file into Illumina's adpc.bin format.

    The adpc.bin format includes::

        x_raw: Raw intensities for allele A
        y_raw: Raw intensities for allele B
        x_norm: Normalized intensities for allele A
        y_norm: Normalized intensities for allele B
        genotype_score: The genotype clustering confidence score
        genotype: The called genotype (0: AA, 1: AB, 2: BB, 3: unknown or missing)
    """

    with AdpcWriter(adpc_file) as fh:
        for record in get_adpc_records(bpm_file, gtc_file):
            fh.write(record)

    if snp_count_file is None:
        snp_count_file = adpc_file.with_suffix(".numSnps.txt")

    with snp_count_file.open("w") as fh:
        num_snps = BeadPoolManifest(bpm_file).num_loci
        fh.write(f"{num_snps}\n")


def get_adpc_records(bpm_file: Path, gtc_file: Path) -> Generator[AdpcRecord, None, None]:
    bpm = BeadPoolManifest(bpm_file)
    gtc = GenotypeCalls(gtc_file)

    x_raw_values = gtc.get_raw_x_intensities()
    y_raw_values = gtc.get_raw_y_intensities()
    x_norm_values, y_norm_values = unzip(gtc.get_normalized_intensities(bpm.normalization_lookups))
    genotype_scores = gtc.get_genotype_scores()
    genotypes = fix_genotype_codes(gtc.get_genotypes())

    for (x_raw, y_raw, x_norm, y_norm, genotype_score, genotype) in zip(
        x_raw_values, y_raw_values, x_norm_values, y_norm_values, genotype_scores, genotypes
    ):
        yield AdpcRecord(x_raw, y_raw, x_norm, y_norm, genotype_score, genotype)


def fix_genotype_codes(genotypes: List[int]) -> Generator[int, None, None]:
    """Convert genotype codes from GTC to Adpc format.

    ===  ====  ===============
    GTC  Adpc    Description
    ===  ====  ===============
    0    3     Unknown/Missing
    1    0     AA
    2    1     AB
    3    2     BB
    ===  ====  ===============
    """
    for genotype in genotypes:
        if genotype in [1, 2, 3]:
            yield genotype - 1
        else:  # if genotype 0 change to 3
            yield 3


if __name__ == "__main__":
    if "snakemake" in locals():
        main(
            Path(snakemake.input.gtc),  # type: ignore # noqa
            Path(snakemake.input.bpm),  # type: ignore # noqa
            Path(snakemake.output.adpc),  # type: ignore # noqa
            Path(snakemake.output.snp_count),  # type: ignore # noqa
        )
    else:
        app()
