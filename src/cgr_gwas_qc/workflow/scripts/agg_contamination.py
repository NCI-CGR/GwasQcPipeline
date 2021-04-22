#!/usr/bin/env python
from pathlib import Path
from typing import List

import pandas as pd
import typer

from cgr_gwas_qc.parsers import plink, verifyidintensity
from cgr_gwas_qc.typing import PathLike

app = typer.Typer(add_completion=False)

DTYPES = {
    "Sample_ID": "string",
    "%Mix": "float",
    "LLK": "float",
    "LLK0": "float",
    "is_contaminated": "boolean",
}


def read(filename: PathLike):
    return pd.read_csv(filename, dtype=DTYPES)


@app.command()
def main(
    contamination_files: List[Path],
    median_intensity_file: Path,
    imiss_file: Path,
    intensity_threshold: int,
    contam_threshold: int,
    outfile: Path,
) -> None:

    df = (
        build(contamination_files, median_intensity_file, imiss_file)
        .pipe(_mask_low_intensity, intensity_threshold)
        .pip(_flag_contaminated, contam_threshold)
    )

    df.reindex(DTYPES.keys()).to_csv(outfile, index=False)


def build(
    contamination_files: List[Path], median_intensity_file: Path, imiss_file: Path
) -> pd.DataFrame:
    return (
        pd.concat(
            [
                _aggregate_contamination(contamination_files).set_index("Sample_ID"),
                pd.read_csv(median_intensity_file).set_index("Sample_ID"),
                plink.read_imiss(imiss_file).rename_axis("Sample_ID"),
            ],
            axis=1,
        )
        .rename_axis("Sample_ID")
        .reset_index()
    )


def _aggregate_contamination(filenames: List[Path]) -> pd.DataFrame:
    return pd.concat([verifyidintensity.read(sample) for sample in filenames], ignore_index=True)


def _mask_low_intensity(df: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Set %Mix to NA if below intensity threshold or not in imiss file"""
    mask = (df.median_intensity < threshold) | df.F_MISS.isna()
    df.loc[mask, "%Mix"] = pd.NA
    return df


def _flag_contaminated(df: pd.DataFrame, threshold: float) -> pd.DataFrame:
    df["is_contaminated"] = False
    mask = df["%Mix"] >= threshold
    df.loc[mask, "is_contaminated"] = True
    return df


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {
            **{
                k: Path(v) if isinstance(v, str) else [Path(f) for f in v]
                for k, v in snakemake.input.items()  # type: ignore # noqa
            },
            **{k: v for k, v in snakemake.params.items()},  # type: ignore # noqa
            "outfile": Path(snakemake.output[0]),  # type: ignore # noqa
        }
        main(**defaults)
    else:
        app()
