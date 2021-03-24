#!/usr/bin/env python
"""
population_qc_table.py
----------------------

Creates a QC summary table for each population. This is done by aggregating major population-level results including:

- Within population subject reslationships (IBS/IBD)
- Principal components analysis
- Autosomal heterozygosity

Output:
    ``population_level/{population}/qc.csv``

    .. csv-table::
        :header: name, description

        population, The population name
        Subject_ID, The Subject's ID
        QC_Family_ID, An arbitrary ID assigned to each related set of subjects.
        relatives, A list of related Subject_IDs concatenated together with a `|`.
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

    - :mod:`cgr_gwas_qc.workflow.scripts.related_subjects`
    - :func:`cgr_gwas_qc.parsers.plink.read_het`
    - :class:`cgr_gwas_qc.parsers.eigensoft.Eigenvec`

"""
import os
from pathlib import Path
from typing import Generator

import pandas as pd
import typer

from cgr_gwas_qc.parsers import eigensoft, plink

app = typer.Typer(add_completion=False)

DTYPES = {
    "population": "string",
    "Subject_ID": "string",
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
}


@app.command()
def main(relatives: Path, pca: Path, autosomal_het: Path, population: str, outfile: Path):
    df = (
        pd.concat([eigensoft.Eigenvec(pca).components, plink.read_het(autosomal_het)], axis=1)
        .rename_axis("Subject_ID")
        .reset_index()
        .assign(population=population)
        .pipe(lambda x: annotate_relations(x, relatives, population))
        # .reindex(DTYPES.keys(), axis=1)
    )

    df.to_csv(outfile, index=False)


def _expand_related(relatives: os.PathLike, population: str) -> Generator[pd.DataFrame, None, None]:
    for fam_id, rels in pd.read_csv(relatives).itertuples(index=False):
        for subject in rels.split("|"):
            yield pd.DataFrame(
                {
                    "Subject_ID": [subject],
                    "QC_Family_ID": [f"{population}_{fam_id}"],
                    "relatives": [rels],
                }
            )
    else:
        yield pd.DataFrame(columns=["Subject_ID", "QC_Family_ID", "relatives"])


def annotate_relations(df: pd.DataFrame, relatives: os.PathLike, population: str) -> pd.DataFrame:
    related_df = pd.concat(_expand_related(relatives, population), ignore_index=True)
    return df.merge(related_df, on="Subject_ID", how="left")


def read_population_qc(filename: os.PathLike) -> pd.DataFrame:
    return pd.read_csv(filename, dtype=DTYPES)


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({"outfile": Path(snakemake.output[0])})  # type: ignore # noqa
        defaults.update({k: v for k, v in snakemake.params.items()})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
