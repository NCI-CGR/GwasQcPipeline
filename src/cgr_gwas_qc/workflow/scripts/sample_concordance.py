#!/usr/bin/env python
"""
Sample Concordance Table
------------------------

``workflow/scripts/sample_concordance.py``

.. csv-table::
    :header: name, dtype, description

    Sample_ID1, string, Sample_ID for the first sample in the pairwise comparison.
    Sample_ID2, string, Sample_ID for the second sample in the pairwise comparison.
    Subject_ID1, string, Subject_ID for the first sample in the pairwise comparison.
    Subject_ID2, string, Subject_ID for the second sample in the pairwise comparison.
    is_internal_control1, string, internal control flag for the first sample in the pairwise comparison.
    is_internal_control2, string, internal control flag for the second sample in the pairwise comparison.
    expected_replicate, boolean, True if the pair of samples are known replicates.
    unexpected_replicate, boolean, True if the pair of samples are from different subjects but look identical.
    PLINK_PI_HAT, float,Proportion IBD i.e. ``P(IBD=2) + 0.5 * P(IBD=1)``
    PLINK_concordance, float, Proportion IBS2 ``IBS2 / (IBS0 + IBS1 + IBS2)``
    PLINK_is_ge_pi_hat, boolean, True if PI_HAT was greater than ``software_params.pi_hat_cutoff``
    PLINK_is_ge_concordance, boolean, True if concordance was greater than ``software_params.dup_concordance_cutoff``
    GRAF_HGMR, float, Homozygous Genotype Mismatch Rate (%)
    GRAF_AGMR, float, All Genotype Mismatch Rate (%)
    GRAF_relationship, string, relationship determined by sample genotypes.
    KING_Kinship, float, Estimated kinship coefficient from the SNP data
    KING_relationship, string, The assigned relationship based on Kinship

References:
    - :mod:`cgr_gwas_qc.workflow.scripts.concordance_table`
    - :mod:`cgr_gwas_qc.parsers.king`
    - :mod:`cgr_gwas_qc.parsers.graf`
"""
from itertools import combinations
from pathlib import Path
from typing import List, Tuple

import pandas as pd
import typer

from cgr_gwas_qc.parsers import graf, king, sample_sheet
from cgr_gwas_qc.typing import PathLike
from cgr_gwas_qc.workflow.scripts import concordance_table

app = typer.Typer(add_completion=False)

DTYPES = {
    "Sample_ID1": "string",
    "Sample_ID2": "string",
    "Subject_ID1": "string",
    "Subject_ID2": "string",
    "is_internal_control1": "boolean",
    "is_internal_control2": "boolean",
    "expected_replicate": "boolean",
    "unexpected_replicate": "boolean",
    "PLINK_PI_HAT": "float",
    "PLINK_concordance": "float",
    "PLINK_is_ge_pi_hat": "boolean",
    "PLINK_is_ge_concordance": "boolean",
    "GRAF_HGMR": "float",
    "GRAF_AGMR": "float",
    "GRAF_relationship": "string",
    "KING_Kinship": "float",
    "KING_relationship": "string",
}


def read(filename: PathLike):
    """Read the sample concordance table

    Returns:
        pd.DataFrame
        - Sample_ID1
        - Sample_ID2
        - Subject_ID1
        - Subject_ID2
        - is_internal_control1
        - is_internal_control2
        - expected_replicate
        - unexpected_replicate
        - PLINK_PI_HAT
        - PLINK_concordance
        - PLINK_is_ge_pi_hat
        - PLINK_is_ge_concordance
        - GRAF_HGMR
        - GRAF_AGMR
        - GRAF_relationship
        - KING_Kinship
        - KING_relationship
    """
    return pd.read_csv(filename, dtype=DTYPES)


@app.command()
def main(
    sample_sheet_csv: Path, plink_file: Path, graf_file: Path, king_file: Path, outfile: Path,
):
    ss = sample_sheet.read(sample_sheet_csv)
    concordance = (
        build(plink_file, graf_file, king_file)
        .pipe(_add_expected_replicates, ss)
        .pipe(_add_unexpected_replicates)
        .pipe(_add_subject, ss)
        .pipe(_add_internal_control, ss)
    )
    concordance.reset_index().reindex(DTYPES.keys(), axis=1).to_csv(outfile, index=False)


def build(plink_file: PathLike, graf_file: PathLike, king_file: PathLike):
    """Build the main concordance table."""
    return (
        _plink(plink_file)
        .join(_graf(graf_file), how="outer")
        .join(_king(king_file), how="outer")
        .rename_axis(["Sample_ID1", "Sample_ID2"])
    )


def _add_subject(df: pd.DataFrame, ss: pd.DataFrame) -> pd.DataFrame:
    """Add Subject IDs for each sample in the pair"""
    s2s = ss.set_index("Sample_ID").Group_By_Subject_ID
    return df.join(s2s.rename_axis("Sample_ID1").rename("Subject_ID1")).join(
        s2s.rename_axis("Sample_ID2").rename("Subject_ID2")
    )


def _add_internal_control(df: pd.DataFrame, ss: pd.DataFrame) -> pd.DataFrame:
    """Add is_internal_control flag for each Sample in pair"""
    s2ic = ss.set_index("Sample_ID").is_internal_control
    return df.join(s2ic.rename_axis("Sample_ID1").rename("is_internal_control1")).join(
        s2ic.rename_axis("Sample_ID2").rename("is_internal_control2")
    )


def _add_expected_replicates(df: pd.DataFrame, ss: pd.DataFrame) -> pd.DataFrame:
    """Flag samples that are expected to be replicates.

    If an expected sample pair is not in the concordance table, then add that
    pair.
    """
    known_replicates = _get_known_replicates(ss)
    df["expected_replicate"] = False
    for pair in known_replicates:
        if pair in df.index:
            df.loc[pair, "expected_replicate"] = True
        else:
            record = pd.Series({"expected_replicate": True}, name=pair)
            df = df.append(record)
    return df


def _add_unexpected_replicates(df: pd.DataFrame) -> pd.DataFrame:
    """Flag pairs of samples that appear to be unexpected replicates.

    Using the different concordance measures, flag a pairs of samples that
    appear to be replicates but are from different subjects.
    """
    df["unexpected_replicate"] = (~df.expected_replicate) & (
        df.PLINK_is_ge_concordance | (df.GRAF_relationship == "ID") | (df.KING_relationship == "ID")
    )
    return df


def _plink(filename: PathLike):
    return (
        concordance_table.read(filename)
        .set_index(["ID1", "ID2"])
        .rename(
            {
                "PI_HAT": "PLINK_PI_HAT",
                "concordance": "PLINK_concordance",
                "is_ge_pi_hat": "PLINK_is_ge_pi_hat",
                "is_ge_concordance": "PLINK_is_ge_concordance",
            },
            axis=1,
        )
    )


def _graf(filename: PathLike):
    return (
        graf.read_relatedness(filename)
        .set_index(["ID1", "ID2"])
        .reindex(["HGMR", "AGMR", "relationship"], axis=1)
        .rename(
            {"HGMR": "GRAF_HGMR", "AGMR": "GRAF_AGMR", "relationship": "GRAF_relationship"}, axis=1,
        )
    )


def _king(filename: PathLike):
    return (
        king.read_kinship(filename)
        .set_index(["ID1", "ID2"])
        .reindex(["Kinship", "relationship"], axis=1)
        .rename({"Kinship": "KING_Kinship", "relationship": "KING_relationship"}, axis=1)
    )


def _get_known_replicates(ss: pd.DataFrame) -> List[Tuple[str, str]]:
    """Create list of pairwise combinations of know replicates.

    This function creates a list of all pairwise combinations of replicate
    IDs. It does this using the following logic::

    - Subset each set of replicates
    - Add all pairwise combinations of replicate IDs to a list

    The ``cgr_sample_sheet.csv`` has a column ``replicate_ids`` which
    contains, for each subject, all of the replicate IDs concatenated together (i.e.,
    "ID1|ID2|ID3"). If there are no replicates than it is ``NA``. When I
    group by ``replicate_ids``, I am essentially subseting subjects that have
    a replicate. I then create all pairwise combinations and add them to a
    list.

    I also make sure to sort ID pairs alphabetically to match the concordance
    table ID order.
    """
    known_replicates = []  # type: ignore
    for _, dd in ss.groupby("replicate_ids"):
        known_replicates.extend(combinations(dd.Sample_ID.sort_values(), 2))
    return known_replicates


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {
            **{k: Path(v) for k, v in snakemake.input.items()},  # type: ignore # noqa
            **{"outfile": Path(snakemake.output[0])},  # type: ignore # noqa
        }
        main(**defaults)
    else:
        app()
