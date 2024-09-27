"""
Aggregate Population Level Concordance
--------------------------------------

.. csv-table::
    :header: name, dtype, description

    population, string, The compared subjects ancestral population.
    Subject_ID1, string, ID1 for the pairwise comparison.
    Subject_ID2, string, ID2 for the pairwise comparison.
    case_control1, string, case/control information for member of the pairwise comparison.
    case_control2, string, case/control information for member of the pairwise comparison.
    PLINK_PI_HAT, float, Proportion IBD i.e. ``P(IBD=2) + 0.5 * P(IBD=1)``
    PLINK_is_ge_pi_hat, boolean, True if PI_HAT was greater than ``software_params.pi_hat_cutoff``
    PLINK_is_ge_concordance, boolean, True if concordance was greater than ``software_params.dup_concordance_cutoff``
"""

from pathlib import Path
from typing import List

import pandas as pd
import typer

from cgr_gwas_qc.parsers import sample_sheet
from cgr_gwas_qc.typing import PathLike
from cgr_gwas_qc.workflow.scripts import concordance_table

app = typer.Typer(add_completion=False)

DTYPES = {
    "population": "string",
    "Subject_ID1": "string",
    "Subject_ID2": "string",
    "case_control1": "string",
    "case_control2": "string",
    "PLINK_PI_HAT": "float",
    "PLINK_is_ge_pi_hat": "boolean",
}


def read(filename: PathLike) -> pd.DataFrame:
    """Reads aggregated population concordance table.

    - population
    - Subject_ID1
    - Subject_ID2
    - case_control1
    - case_control2
    - PLINK_PI_HAT
    - PLINK_is_ge_pi_hat
    """
    return pd.read_csv(filename, dtype=DTYPES)


@app.command()
def main(sample_sheet_csv: Path, ibd_files: List[Path], outfile: Path):
    if not ibd_files:
        pd.DataFrame(columns=DTYPES.keys()).to_csv(outfile, index=False)
        return None

    ss = (
        sample_sheet.read(sample_sheet_csv)
        .reindex(["Group_By_Subject_ID", "case_control"], axis=1)
        .drop_duplicates()
        .set_index("Group_By_Subject_ID")
        .squeeze()
    )

    df = (
        pd.concat([_read_concordance_table(filename) for filename in ibd_files], ignore_index=True)
        .merge(ss.rename_axis("Subject_ID1").rename("case_control1"), on="Subject_ID1", how="left")
        .merge(ss.rename_axis("Subject_ID2").rename("case_control2"), on="Subject_ID2", how="left")
    )

    df.reindex(DTYPES.keys(), axis=1).to_csv(outfile, index=False)


def _read_concordance_table(filename: Path):
    population = filename.parent.stem
    return (
        concordance_table.read(filename)
        .rename(
            {
                "ID1": "Subject_ID1",
                "ID2": "Subject_ID2",
                "PI_HAT": "PLINK_PI_HAT",
                "is_ge_pi_hat": "PLINK_is_ge_pi_hat",
            },
            axis=1,
        )
        .assign(population=population)
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        ibd_files = snakemake.input.ibd_files  # type: ignore # noqa
        if isinstance(ibd_files, str):
            filenames = [Path(ibd_files)]
        else:
            filenames = [Path(x) for x in ibd_files]

        defaults = {
            "sample_sheet_csv": Path(snakemake.input.sample_sheet_csv),  # type: ignore # noqa
            "ibd_files": filenames,  # type: ignore # noqa
            "outfile": Path(snakemake.output[0]),  # type: ignore # noqa
        }
        main(**defaults)
    else:
        app()
