"""Population level summary table with the following columns:

.. csv-table::
    :header: name, dtype, description

    "Subject_ID", str, Subject ID
    "Sample_ID", str, Sample ID
    "case_control", str, Case/Control Status
    "population", str, Population name
    "PC1", float, Principal component 1
    "PC2", float, Principal component 2
    "PC3", float, Principal component 3
    "PC4", float, Principal component 4
    "PC5", float, Principal component 5
    "PC6", float, Principal component 6
    "PC7", float, Principal component 7
    "PC8", float, Principal component 8
    "PC9", float, Principal component 9
    "PC10", float, Principal component 10
    "O(HOM)", int, Observed number of homozygotes
    "E(HOM)", int, Expected number of homozygotes
    "N(NM)", int, Number of non-missing autosomal genotypes
    "F", float, Method-of-moments F coefficient estimate

References:

    - :func:`cgr_gwas_qc.parsers.plink.read_het`
    - :class:`cgr_gwas_qc.parsers.eigensoft.Eigenvec`

"""
from dataclasses import dataclass
from itertools import chain
from pathlib import Path
from typing import Generator

import pandas as pd
import typer

from cgr_gwas_qc.parsers import eigensoft, plink

app = typer.Typer(add_completion=False)

col_order = [
    "Subject_ID",
    "Sample_ID",
    "case_control",
    "population",
    "PC1",
    "PC2",
    "PC3",
    "PC4",
    "PC5",
    "PC6",
    "PC7",
    "PC8",
    "PC9",
    "PC10",
    "O(HOM)",
    "E(HOM)",
    "N(NM)",
    "F",
]


@app.command()
def main(sample_qc: Path, results: Path, controls: Path, outfile: Path):
    df = (
        build_table(results, controls)
        .pipe(add_metadata, filename=sample_qc)
        .reindex(col_order, axis=1)
    )

    df.to_csv(outfile, index=False)


def build_table(results: Path, controls: Path) -> pd.DataFrame:
    read_funcs = {
        "eigenvec": _eigenvec,
        "het": _het,
    }

    data = [
        read_funcs[popfile.suffix](popfile)
        for popfile in chain(extract_files(results), extract_files(controls))
        if popfile.suffix in read_funcs
    ]

    if not data:
        return pd.DataFrame()

    return (
        pd.concat(data)
        .set_index(["Subject_ID", "population"])
        .stack()  # this gets ride of the NA and merges records together for each subject*population.
        .unstack()
        .reset_index()
        .sort_values(["population", "Subject_ID"])
        .astype({"O(HOM)": int, "E(HOM)": int, "N(NM)": int})
    )


def add_metadata(df: pd.DataFrame, filename: Path):
    metadata = (
        pd.read_csv(filename)
        .query("is_subject_representative")
        .reindex(["Group_By_Subject_ID", "Sample_ID", "case_control"], axis=1)
        .rename({"Group_By_Subject_ID": "Subject_ID"}, axis=1)
    )
    try:
        return df.merge(metadata, on="Subject_ID", how="left")
    except KeyError:
        return pd.DataFrame()


@dataclass
class PopFile:
    population: str
    path: Path
    suffix: str


def extract_files(filelist: Path) -> Generator[PopFile, None, None]:
    for filename in filelist.read_text().strip().splitlines():
        path = Path(filename)
        suffix = path.suffix.lstrip(".")
        population = path.parent.name
        yield PopFile(population, path, suffix)


def _eigenvec(popfile: PopFile) -> pd.DataFrame:
    return (
        eigensoft.Eigenvec(popfile.path)
        .components.reset_index()
        .rename({"ID": "Subject_ID"}, axis=1)
        .assign(population=popfile.population)
    )


def _het(popfile: PopFile) -> pd.DataFrame:
    return (
        plink.read_het(popfile.path)
        .reset_index()
        .rename({"ID": "Subject_ID"}, axis=1)
        .assign(population=popfile.population)
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({"outfile": Path(snakemake.output[0])})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
