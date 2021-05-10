"""
agg_population_qc_tables.py
---------------------------

Aggregate per population summary tables into a single table. Then add extra
metadata from the sample qc table.

Output:
    ``population_level/population_qc.csv``

    .. csv-table::
        :header: name, description

        population, The population name
        Subject_ID, The Subject's ID
        Sample_ID, Sample ID
        case_control, Case/Control Status
        QC_Family_ID, An arbitrary ID assigned to each related set of subjects.
        relatives, A list of related Subject_IDs concatenated together with a '|'.
        PC1, Principal component 1
        PC2, Principal component 2
        PC3, Principal component 3
        PC4, Principal component 4
        PC5, Principal component 5
        PC6, Principal component 6
        PC7, Principal component 7
        PC8, Principal component 8
        PC9, Principal component 9
        PC10, Principal component 10
        O_HOM, Observed number of homozygotes
        E_HOM, Expected number of homozygotes
        N_NM, Number of non-missing autosomal genotypes
        F, Method-of-moments F coefficient estimate

References:

    - :mod:`cgr_gwas_qc.workflow.scripts.sample_qc_table`
    - :mod:`cgr_gwas_qc.workflow.scripts.related_subjects`
    - :func:`cgr_gwas_qc.parsers.plink.read_het`
    - :class:`cgr_gwas_qc.parsers.eigensoft.Eigenvec`

"""
from pathlib import Path
from typing import List

import pandas as pd
import typer

from cgr_gwas_qc.reporting import CASE_CONTROL_DTYPE
from cgr_gwas_qc.workflow.scripts import population_qc_table

app = typer.Typer(add_completion=False)

DTYPES = {
    "population": "string",
    "Subject_ID": "string",
    "Sample_ID": "string",
    "case_control": CASE_CONTROL_DTYPE,
    "QC_Family_ID": "string",
    "relatives": "string",
    "PC1": "float",
    "PC2": "float",
    "PC3": "float",
    "PC4": "float",
    "PC5": "float",
    "PC6": "float",
    "PC7": "float",
    "PC8": "float",
    "PC9": "float",
    "PC10": "float",
    "O_HOM": "UInt32",
    "E_HOM": "UInt32",
    "N_NM": "UInt32",
    "F": "float",
    "is_extreme_autosomal_heterozygosity": "boolean",
}


@app.command()
def main(subject_qc_table: Path, population_qc_tables: List[Path], outfile: Path):
    if not population_qc_tables:
        # No populations to aggregate save an empty table
        pd.DataFrame(columns=DTYPES.keys()).to_csv(outfile, index=False)
        return None

    df = (
        aggregate_qc_tables(population_qc_tables)
        .pipe(add_metadata, filename=subject_qc_table)
        .reindex(DTYPES.keys(), axis=1)
    )

    df.to_csv(outfile, index=False)


def aggregate_qc_tables(population_files: List[Path]) -> pd.DataFrame:
    return pd.concat(
        [population_qc_table.read(Path(filename)) for filename in population_files],
        ignore_index=True,
    )


def add_metadata(df: pd.DataFrame, filename: Path):
    metadata = (
        pd.read_csv(filename)
        .reindex(["Group_By_Subject_ID", "Sample_ID", "case_control"], axis=1)
        .rename({"Group_By_Subject_ID": "Subject_ID"}, axis=1)
    )
    try:
        return df.merge(metadata, on="Subject_ID", how="left")
    except KeyError:
        return pd.DataFrame()


def read_agg_population_qc_tables(filename: Path):
    return pd.read_csv(filename, dtype=DTYPES)


if __name__ == "__main__":
    if "snakemake" in locals():
        qc_files = snakemake.input.population_qc_tables  # type: ignore # noqa
        if isinstance(qc_files, str):
            filenames = [Path(qc_files)]
        else:
            filenames = [Path(x) for x in qc_files]

        defaults = {
            "subject_qc_table": Path(snakemake.input.subject_qc_table),  # type: ignore # noqa
            "population_qc_tables": filenames,
            "outfile": Path(snakemake.output[0]),  # type: ignore # noqa
        }
        main(**defaults)
    else:
        app()
