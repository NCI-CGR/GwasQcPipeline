"""Creates a SNP QC Summary Table.

Columns:
    - array_snp_id (str)
    - thousand_genomes_snp_id (str)
    - chromosome (int, str)
    - snp_cr1 (float)
    - snp_cr1_removed (bool)
    - snp_cr2 (float)
    - snp_cr2_removed (bool)
"""
from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False)


@app.command()
def main(
    initial: Path = typer.Option(..., help="Initial lmiss file."),
    cr1: Path = typer.Option(..., help="lmiss file after CR1 filters."),
    cr2: Path = typer.Option(..., help="lmiss file after CR2 filters."),
    thousand_genomes: Path = typer.Option(..., help="Mapping of Array IDs to thousand genome IDs."),
    outfile: Path = typer.Option(..., help="Output CSV."),
):

    build_table(initial, cr1, cr2, thousand_genomes).to_csv(outfile)


def build_table(initial, cr1, cr2, thousand_genomes):
    return (
        _read_lmiss(initial, "call_rate_initial")
        .merge(_read_lmiss(cr1, "cr1"), how="outer")
        .merge(_read_lmiss(cr2, "cr2"), how="outer")
        .merge(_read_thousand_genomes(thousand_genomes), how="outer")
        .fillna({"snp_cr1_removed": True, "snp_cr2_removed": True})
        .reindex(
            [
                "array_snp_id",
                "thousand_genomes_snp_id",
                "chromosome",
                "snp_cr1",
                "snp_cr1_removed",
                "snp_cr2",
                "snp_cr2_removed",
            ],
            axis=1,
        )
    )


def _read_lmiss(filename: Path, prefix: str) -> pd.DataFrame:
    df = (
        pd.read_csv(filename, delim_whitespace=True)
        .assign(**{f"snp_{prefix}": lambda x: 1 - x.F_MISS})
        .rename({"SNP": "array_snp_id", "CHR": "chromosome"}, axis=1)
        .reindex(["array_snp_id", "chromosome", f"snp_{prefix}"], axis=1)
    )

    if prefix != "call_rate_initial":
        df[f"snp_{prefix}_removed"] = False

    return df


def _read_thousand_genomes(filename: Path):
    return pd.read_csv(filename).rename(
        {"array_id": "array_snp_id", "thousand_genomes_id": "thousand_genomes_snp_id"}, axis=1
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type:ignore # noqa
        defaults.update({"outfile": Path(snakemake.output[0])})  # type:ignore # noqa
        main(**defaults)
    else:
        app()
