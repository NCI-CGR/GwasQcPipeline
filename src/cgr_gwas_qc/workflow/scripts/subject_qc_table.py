#!/usr/bin/env python
"""
Internal Subject QC Report
--------------------------

``workflow/scripts/subject_qc_table.py``

.. csv-table::
    :header: name, dtype, description

    Group_By_Subject_ID, string, The subject identifier used by the workflow
    Sample_ID, string, The sample identifier used by the workflow.
    subject_dropped_from_study, boolean, True if there are no samples passing QC for this subject.
    case_control, CASE_CONTROL_DTYPE, Phenotype status {Case, Control, QC, Unknown}
    is_unexpected_replicate, boolean, True if two subjects had high concordance.
    unexpected_replicate_ids, string, Concatenated Sample_IDs that are unexpected replicates.
    expected_sex, SEX_DTYPE, The expected sex from the provided sample sheet.
    predicted_sex, SEX_DTYPE, The predicted sex based on X chromosome heterozygosity.
    X_inbreeding_coefficient, float, The X chromosome F coefficient.
    is_sex_discordant, boolean, True if expected and predicted sex did not match.
    AFR, float, Proportion African ancestry.
    EUR, float, Proportion European ancestry.
    ASN, float, Proportion East Asian ancestry.
    Ancestry, category, Assigned ancestral population (GRAF).
"""
from itertools import product
from pathlib import Path
from typing import Dict, Iterable, Tuple

import networkx as nx
import pandas as pd
import typer

from cgr_gwas_qc.reporting import CASE_CONTROL_DTYPE, SEX_DTYPE
from cgr_gwas_qc.typing import PathLike
from cgr_gwas_qc.workflow.scripts import sample_concordance, sample_qc_table

app = typer.Typer(add_completion=False)

DTYPES = {  # Header for main QC table
    "Group_By_Subject_ID": "string",
    "Sample_ID": "string",
    "case_control": CASE_CONTROL_DTYPE,
    "is_unexpected_replicate": "boolean",
    "unexpected_replicate_ids": "string",
    "expected_sex": SEX_DTYPE,
    "predicted_sex": SEX_DTYPE,
    "X_inbreeding_coefficient": "float",
    "is_sex_discordant": "boolean",
    "AFR": "float",
    "EUR": "float",
    "ASN": "float",
    "Ancestry": "category",
}


def read(filename: PathLike) -> pd.DataFrame:
    """Read the subject level QC table.

    Returns:
        pd.DataFrame:
        - Group_By_Subject_ID
        - Sample_ID
        - case_control
        - is_unexpected_replicate
        - unexpected_replicate_ids
        - expected_sex
        - predicted_sex
        - X_inbreeding_coefficient
        - is_sex_discordant
        - AFR
        - EUR
        - ASN
        - Ancestry
    """
    return pd.read_csv(filename, dtype=DTYPES)


@app.command()
def main(
    sample_qc_csv: Path, sample_concordance_csv: Path, outfile: Path,
):
    (
        sample_qc_table.read(sample_qc_csv)
        .pipe(_fix_hyphen_in_ancestry_name)
        .pipe(_sample_qc_to_subject_qc)
        .pipe(_add_unexpected_replicate_ids, sample_concordance_csv)
        .reindex(DTYPES.keys(), axis=1)
        .to_csv(outfile, index=False)
    )


#Needed with graf-pop implementarion that returns Asian-Pacific_Islander which the hyphen isn't supported by snakemake wildcares
def _fix_hyphen_in_ancestry_name(df: pd.DataFrame) -> pd.DataFrame:
    df['Ancestry'] = df['Ancestry'].str.replace('-', '_')
    return (
        df
    )


def _sample_qc_to_subject_qc(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.query("is_subject_representative")
        .reindex(DTYPES.keys(), axis=1)
        .dropna(how="all", axis=1)
    )


def _add_unexpected_replicate_ids(df: pd.DataFrame, sample_concordance_csv: Path) -> pd.DataFrame:
    """Adds information about unexpected replicates.

    Adds ``is_unexpected_replicate`` and ``unexpected_replicate_ids``.
    """
    concordance = sample_concordance.read(sample_concordance_csv)
    flags = (
        concordance.melt(
            id_vars=["is_unexpected_replicate"],
            value_vars=["Subject_ID1", "Subject_ID2"],
            var_name="To_Drop",
            value_name="Group_By_Subject_ID",
        )
        .drop("To_Drop", axis=1)
        .groupby("Group_By_Subject_ID")
        .max()  # Flag a sample as True if it is True for any comparison.
        .astype("boolean")
    )

    unexpected_ids = (
        _connected_ids(
            concordance.set_index(["Subject_ID1", "Subject_ID2"])
            .query("is_unexpected_replicate")
            .index.values
        )
        .rename_axis("Subject_ID")
        .rename("unexpected_replicate_ids")
        .astype("string")
    )

    return df.merge(flags.join(unexpected_ids), on="Group_By_Subject_ID", how="left").fillna(
        {"is_unexpected_replicate": False}
    )


def _connected_ids(ids: Iterable[Tuple[str, str]]) -> pd.Series:
    """Create groups of connected IDs.

    Given an iterable of tuples, where each tuple represents an edge in a graph. Builds
    the full graph and then for each subgraph creates a list of connected ids
    in the form `ID1|ID2|...`.

    Example
    -------
    >>> dict(_connected_ids([("one", "two"), ("two", "three"), ("four", "five")]))
    ... {
            "one": "one|two|three",
            "two": "one|two|three",
            "three": "one|two|three",
            "four": "four|five",
            "five": "four|five",
        }
    """
    G = nx.Graph()
    G.add_edges_from(ids)
    groups: Dict[str, str] = {}

    for subgraph in sorted(nx.connected_components(G), key=len):
        ids_string = "|".join(sorted(subgraph))
        groups.update(product(subgraph, [ids_string]))

    return pd.Series(groups, dtype="string")


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {
            **{k: Path(v) for k, v in snakemake.input.items()},  # type: ignore # noqa
            **{k: v for k, v in snakemake.params.items()},  # type: ignore # noqa
            "outfile": Path(snakemake.output[0]),  # type: ignore # noqa
        }
        main(**defaults)
    else:
        app()
