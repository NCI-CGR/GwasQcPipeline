"""Summarize IBS/IBD concordance.

Reads the plink `.genome` format file and calculates concordance as IBS2 /
(IBS0 + IBS1 + IBS2).

.. csv-table::
    :header: name, dtype, description

    ID1, string, ID1 for the pairwise comparison.
    ID2, string, ID2 for the pairwise comparison.
    PI_HAT, float, Proportion IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)
    concordance, float, Proportion IBS2 ``IBS2 / (IBS0 + IBS1 + IBS2)``
    is_ge_pi_hat, boolean, True if PI_HAT was greater than `software_params.pi_hat_cutoff`
    is_ge_concordance, boolean, True if concordance was greater than `software_params.dup_concordance_cutoff`

"""
import os
from pathlib import Path

import pandas as pd
import typer

from cgr_gwas_qc.parsers import plink

app = typer.Typer(add_completion=False)

DTYPES = {
    "ID1": "string",
    "ID2": "string",
    "PI_HAT": "float",
    "concordance": "float",
    "is_ge_pi_hat": "boolean",
    "is_ge_concordance": "boolean",
}


@app.command()
def main(filename: Path, concordance_threshold: float, pi_hat_threshold: float, outfile: Path):
    create_concordance_table(filename, concordance_threshold, pi_hat_threshold).to_csv(outfile)


def read_concordance(filename: os.PathLike) -> pd.DataFrame:
    return pd.read_csv(filename, dtype=DTYPES)


def create_concordance_table(
    filename: Path, concordance_threshold: float, pi_hat_threshold: float
) -> pd.DataFrame:
    return (
        plink.read_genome(filename)
        .assign(is_ge_pi_hat=lambda x: x.PI_HAT >= pi_hat_threshold)
        .assign(concordance=lambda x: x.IBS2 / (x.IBS0 + x.IBS1 + x.IBS2))
        .assign(is_ge_concordance=lambda x: x.concordance >= concordance_threshold)
        .reindex(DTYPES.keys(), axis=1)
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({"filename": Path(snakemake.input[0])})  # type: ignore # noqa
        defaults.update({k: float(v) for k, v in snakemake.params.items()})  # type: ignore # noqa
        defaults.update({"outfile": Path(snakemake.output[0])})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
