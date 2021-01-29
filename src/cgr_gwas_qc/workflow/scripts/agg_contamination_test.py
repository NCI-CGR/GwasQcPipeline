#!/usr/bin/env python
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import typer

from cgr_gwas_qc.parsers.verifyidintensity import parse_sample_contamination_output

app = typer.Typer(add_completion=False)


@app.command()
def main(
    contamination_files: List[Path] = typer.Argument(
        ..., help="List of Paths with sample contamination estimates."
    ),
    median_intensity_file: Path = typer.Argument(
        ..., help="Path to aggregated median IDAT intensities.", exists=True, readable=True
    ),
    imiss_file: Path = typer.Argument(
        ...,
        help="Path to sample based missing report from lvl2 call rate.",
        exists=True,
        readable=True,
    ),
    intensity_threshold: int = typer.Argument(..., help="IDAT intensity threshold."),
    out_file: Path = typer.Argument(
        ..., help="Path to output file.", file_okay=True, writable=True
    ),
) -> None:
    aggregated_sample_contamination = pd.concat(
        [parse_sample_contamination_output(sample) for sample in contamination_files],
        ignore_index=True,
    )
    median_intensity = pd.read_csv(median_intensity_file)
    aggregated_sample_missingness = parse_imiss(imiss_file)

    df = aggregated_sample_contamination.merge(median_intensity, on="Sample_ID").merge(
        aggregated_sample_missingness, on="Sample_ID"
    )

    # Set %Mix to NA if below intensity threshold or not in imiss file
    mask = (df.median_intensity < intensity_threshold) & df.F_MISS.isnull()
    df.loc[mask, "%Mix"] = np.nan

    # Only output a subset of columns. Mostly only care about ``Sample_ID`` and ``%Mix``.
    df[["Sample_ID", "%Mix", "LLK", "LLK0"]].to_csv(out_file, index=False)


def parse_imiss(file_name: Path):
    return pd.read_csv(file_name, delim_whitespace=True).rename({"IID": "Sample_ID"}, axis=1)


if __name__ == "__main__":
    if "snakemake" in locals():
        main(
            [Path(x) for x in snakemake.input.contamination],  # type: ignore # noqa
            Path(snakemake.input.median_idat_intensity),  # type: ignore # noqa
            Path(snakemake.input.imiss),  # type: ignore # noqa
            snakemake.params.intensity_threshold,  # type: ignore # noqa
            Path(snakemake.output[0]),  # type: ignore # noqa
        )
    else:
        app()
