"""
Internal Population QC Report
-----------------------------

**Script**: ``agg_population_qc_tables.py``

Aggregate per population summary tables into a single table. Then add extra
metadata from the sample qc table.

**Output**: ``population_level/population_qc.csv``

    .. csv-table::
        :header: name, dtype, description

        population, string, The population name
        Subject_ID, string, The subject identifier used by the workflow
        Sample_ID, string, The sample identifier used by the workflow.
        case_control, CASE_CONTROL_DTYPE, Phenotype status [``Case`` | ``Control`` | ``QC`` | ``Unknown``].
        QC_Family_ID, string, An arbitrary ID assigned to each related set of subjects.
        relatives, string, A list of related Subject_IDs concatenated together with a '|'.
        PC1, float, Principal component 1
        PC2, float, Principal component 2
        PC3, float, Principal component 3
        PC4, float, Principal component 4
        PC5, float, Principal component 5
        PC6, float, Principal component 6
        PC7, float, Principal component 7
        PC8, float, Principal component 8
        PC9, float, Principal component 9
        PC10, float, Principal component 10
        O_HOM, int, Observed number of homozygotes
        E_HOM, int, Expected number of homozygotes
        N_NM, int, Number of non-missing autosomal genotypes
        F, float, Method-of-moments F coefficient estimate

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
from cgr_gwas_qc.workflow.scripts import population_qc_table, subject_qc_table

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
    # i217####
    gwas_df = create_gwas(df)
    gwas_outfile = "subject_level/gwas.txt"
    gwas_df.to_csv(gwas_outfile, sep=" ", index=False)
    print(gwas_outfile + " complete")
    # i217####


def aggregate_qc_tables(population_files: List[Path]) -> pd.DataFrame:
    return pd.concat(
        [population_qc_table.read(Path(filename)) for filename in population_files],
        ignore_index=True,
    )


def add_metadata(df: pd.DataFrame, filename: Path):
    metadata = (
        subject_qc_table.read(filename)
        .reindex(["Group_By_Subject_ID", "Sample_ID", "case_control"], axis=1)
        .rename({"Group_By_Subject_ID": "Subject_ID"}, axis=1)
    )
    try:
        return df.merge(metadata, on="Subject_ID", how="left")
    except KeyError:
        return pd.DataFrame()


def read_agg_population_qc_tables(filename: Path):
    return pd.read_csv(filename, dtype=DTYPES)


# i217####
def create_gwas(df1):
    temp_df = df1[["Sample_ID", "case_control"]].copy()
    temp_df["case_control"] = temp_df["case_control"].str.lower()
    temp_df = temp_df.loc[temp_df.case_control.str.contains("case|control"), :]
    temp_df["case"] = temp_df["case_control"].map({"case": "2", "control": "1"})
    temp_df["IID"] = temp_df["Sample_ID"]
    temp_df.rename(columns={"Sample_ID": "FID"}, inplace=True)
    temp_df = temp_df[["FID", "IID", "case", "case_control"]]
    gwas_df = temp_df.dropna().copy()
    return gwas_df


# i217####

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
