"""
python gtc2plink.py /path/to/file.gtc /path/to/manifest.bpm /path/to/out.adpc.bin
also outputs a txt file with one line that has the number of markers:
/path/to/out.adpc.bin.numSnps.txt

written by Eric Karlins with help from Bin Zhu

note that I'm using a version of IlluminaBeadArrayFiles.py that I edited to address ZeroDivisionError in the function get_normalized_intensities
"""

from pathlib import Path
from typing import Generator, List, Optional

import typer
from more_itertools import unzip

from cgr_gwas_qc.parsers.adpc import AdpcRecord, AdpcWriter
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
    adpc_file: Path = typer.Argument(..., help="Output adpc file.", file_okay=True, writable=True),
    snp_count_file: Optional[Path] = typer.Argument(
        None, help="Save the number of snps.", file_okay=True, writable=True
    ),
) -> None:

    with AdpcWriter(adpc_file) as fh:
        for record in get_adpc_records(bpm_file, gtc_file):
            fh.write(record)

    if snp_count_file is None:
        snp_count_file = adpc_file.with_suffix(".numSnps.txt")

    with snp_count_file.open("w") as fh:
        num_snps = BeadPoolManifest(bpm_file).num_loci
        fh.write(f"{num_snps}\n")


def get_adpc_records(bpm_file: Path, gtc_file: Path) -> Generator["AdpcRecord", None, None]:
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
    app()
