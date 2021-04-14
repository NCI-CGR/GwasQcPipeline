"""
Aggregate Population Level Concordance
--------------------------------------

.. csv-table::
    :header: name, dtype, description

    population, string, The compared subjects ancestral population.
    Subject_ID1, string, ID1 for the pairwise comparison.
    Subject_ID2, string, ID2 for the pairwise comparison.
    PI_HAT, float, Proportion IBD i.e. ``P(IBD=2) + 0.5 * P(IBD=1)``
    concordance, float, Proportion IBS2 ``IBS2 / (IBS0 + IBS1 + IBS2)``
    is_ge_pi_hat, boolean, True if PI_HAT was greater than ``software_params.pi_hat_cutoff``
    is_ge_concordance, boolean, True if concordance was greater than ``software_params.dup_concordance_cutoff``
"""
from pathlib import Path
from typing import List

import pandas as pd
import typer

from cgr_gwas_qc.typing import PathLike
from cgr_gwas_qc.workflow.scripts import concordance_table

app = typer.Typer(add_completion=False)

DTYPES = {
    "population": "string",
    "Subject_ID1": "string",
    "Subject_ID2": "string",
    "PI_HAT": "float",
    "concordance": "float",
    "is_ge_pi_hat": "boolean",
    "is_ge_concordance": "boolean",
}


def read(filename: PathLike) -> pd.DataFrame:
    """Reads aggregated population concordance table.

    - ``population``
    - ``Subject_ID1``
    - ``Subject_ID2``
    - ``PI_HAT``
    - ``concordance``
    - ``is_ge_pi_hat``
    - ``is_ge_concordance``
    """
    return pd.read_csv(filename, dtype=DTYPES)


@app.command()
def main(filenames: List[Path], outfile: Path):
    if filenames:
        df = pd.concat(
            [_read_concordance_table(filename) for filename in filenames], ignore_index=True
        )
    else:
        df = pd.DataFrame([], columns=DTYPES.keys(), dtype=DTYPES)

    df.reindex(DTYPES.keys(), axis=1).to_csv(outfile, index=False)


def _read_concordance_table(filename: Path):
    population = filename.parent.stem
    return (
        concordance_table.read(filename)
        .rename({"ID1": "Subject_ID1", "ID2": "Subject_ID2"}, axis=1)
        .assign(population=population)
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {
            "filenames": [Path(filename) for filename in snakemake.input],  # type: ignore # noqa
            "outfile": Path(snakemake.output[0]),  # type: ignore # noqa
        }
        main(**defaults)
    else:
        app()
