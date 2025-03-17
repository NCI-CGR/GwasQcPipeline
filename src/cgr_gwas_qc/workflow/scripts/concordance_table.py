"""
Base Concordance Table
----------------------

Reads the plink ``.genome`` format file and calculates concordance as
``IBS2 / (IBS0 + IBS1 + IBS2)``.

.. csv-table::
    :header: name, dtype, description

    ID1, string, ID1 for the pairwise comparison.
    ID2, string, ID2 for the pairwise comparison.
    PI_HAT, float, Proportion IBD i.e. ``P(IBD=2) + 0.5 * P(IBD=1)``
    concordance, float, Proportion IBS2 ``IBS2 / (IBS0 + IBS1 + IBS2)``
    is_ge_pi_hat, boolean, True if PI_HAT was greater than ``software_params.pi_hat_cutoff``
    is_ge_concordance, boolean, True if concordance was greater than ``software_params.dup_concordance_cutoff``

"""

import concurrent.futures
import functools
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import typer

from cgr_gwas_qc.parsers import plink
from cgr_gwas_qc.typing import PathLike

app = typer.Typer(add_completion=False)

DTYPES = {
    "ID1": "category",
    "ID2": "category",
    "PI_HAT": "float",
    "concordance": "float",
    "is_ge_pi_hat": "boolean",
    "is_ge_concordance": "boolean",
}


@app.command()
def main(
    filename: Path,
    concordance_threshold: float,
    pi_hat_threshold: float,
    outfile: Path,
    threads: int = 1,
):
    build(filename, concordance_threshold, pi_hat_threshold, outfile, threads)


def read(filename: PathLike) -> pd.DataFrame:
    """Reads concordance table.

    Returns:
        A table with:

        - ``ID1``
        - ``ID2``
        - ``PI_HAT``
        - ``concordance``
        - ``is_ge_pi_hat``
        - ``is_ge_concordance``
    """
    return pd.read_csv(filename, dtype=DTYPES)


def _sort_ids(x: pd.DataFrame):
    """Sort IDs alphanumerically."""
    x.IID1, x.IID2 = np.where(x.IID1 < x.IID2, [x.IID1, x.IID2], [x.IID2, x.IID1])
    x.rename(columns={"IID1": "ID1", "IID2": "ID2"}, inplace=True)
    return x


def prep_concordance_table(pi_hat_threshold: float, concordance_threshold: float, x: pd.DataFrame):
    """Prepares the concordance_table by appending new columns based on thresholds."""
    return (
        _sort_ids(x)
        .assign(is_ge_pi_hat=lambda x: x.PI_HAT >= pi_hat_threshold)
        .assign(concordance=lambda x: x.IBS2 / (x.IBS0 + x.IBS1 + x.IBS2))
        .assign(is_ge_concordance=lambda x: x.concordance >= concordance_threshold)
        .reindex(DTYPES.keys(), axis=1)
    )


def process_genome_chunk(pi_hat_threshold: float, concordance_threshold: float, chunk):
    """Processes each chunk of genome file."""
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    prep_concordance_table(pi_hat_threshold, concordance_threshold, chunk).to_csv(
        temp_file.name, index=False, header=False
    )
    return temp_file.name


def build(
    filename: Path,
    concordance_threshold: float,
    pi_hat_threshold: float,
    outfile: Path,
    threads: int,
):

    plink_genome_file = plink.read_genome(
        filename, required_cols=["IID1", "IID2", "PI_HAT", "IBS0", "IBS1", "IBS2"]
    )

    temp_files = []

    header = tempfile.NamedTemporaryFile(delete=False)
    temp_files.append(header.name)
    pd.DataFrame(columns=DTYPES.keys()).to_csv(header.name, index=False)

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        results = executor.map(
            functools.partial(process_genome_chunk, pi_hat_threshold, concordance_threshold),
            plink_genome_file,
        )

    temp_files = temp_files + list(results)
    with open(outfile, "wb") as outputfile:
        subprocess.run(["cat"] + temp_files, stdout=outputfile)

    for f in temp_files:
        Path(f).unlink()


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({"filename": Path(snakemake.input[0])})  # type: ignore # noqa
        defaults.update({k: float(v) for k, v in snakemake.params.items()})  # type: ignore # noqa
        defaults.update({"outfile": Path(snakemake.output[0])})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
