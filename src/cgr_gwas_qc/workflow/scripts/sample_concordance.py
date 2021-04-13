#!/usr/bin/env python
"""
Sample Concordance Tables
-------------------------

This script uses sample concordance based on pairwise IBD estimations. We
expect samples from the same subject to show high concordance while samples
from different subjects should not be related. This script outputs 4 tables:

Known Replicates
++++++++++++++++

- ``sample_level/concordance/KnownReplicates.csv`` (Full Table)
- ``sample_level/concordance/InternalQcKnown.csv`` (QC Samples Only)
- ``sample_level/concordance/StudySampleKnown.csv`` (Study Samples Only)

These tables have the following format:

.. csv-table::
    :header: name, dtype, description

    Subject_ID, string, The Subject_ID for the replicate
    Sample_ID1, string, Sample_ID1 for the pairwise comparison.
    Sample_ID2, string, Sample_ID2 for the pairwise comparison.
    PI_HAT, float, Proportion IBD i.e. ``P(IBD=2) + 0.5*P(IBD=1)``
    concordance, float, Proportion IBS2 ``IBS2 / (IBS0 + IBS1 + IBS2)``
    is_ge_pi_hat, boolean, True if PI_HAT was greater than ``software_params.pi_hat_cutoff``
    is_ge_concordance, boolean, True if concordance was greater than ``software_params.dup_concordance_cutoff``

Unknown Replicates
++++++++++++++++++

- ``sample_level/concordance/UnknownReplicates.csv`` (Full Table)

This table has the following format:

.. csv-table::
    :header: name, dtype, description

    Subject_ID1, string, The Subject_ID for ``Sample_ID1``
    Subject_ID2, string, The Subject_ID for ``Sample_ID2``
    Sample_ID1, string, Sample_ID1 for the pairwise comparison.
    Sample_ID2, string, Sample_ID2 for the pairwise comparison.
    PI_HAT, float, Proportion IBD i.e. ``P(IBD=2) + 0.5*P(IBD=1)``
    concordance, float, Proportion IBS2 ``IBS2 / (IBS0 + IBS1 + IBS2)``
    is_ge_pi_hat, boolean, True if PI_HAT was greater than ``software_params.pi_hat_cutoff``
    is_ge_concordance, boolean, True if concordance was greater than ``software_params.dup_concordance_cutoff``

"""
from itertools import combinations
from pathlib import Path
from typing import List, Tuple

import pandas as pd
import typer

from cgr_gwas_qc.parsers import sample_sheet
from cgr_gwas_qc.typing import PathLike
from cgr_gwas_qc.workflow.scripts import concordance_table

app = typer.Typer(add_completion=False)

KNOWN_DTYPES = {
    "Subject_ID": "string",
    "Sample_ID1": "string",
    "Sample_ID2": "string",
    "PI_HAT": "float",
    "concordance": "float",
    "is_ge_pi_hat": "boolean",
    "is_ge_concordance": "boolean",
}

UNKNOWN_DTYPES = {
    "Subject_ID1": "string",
    "Subject_ID2": "string",
    "Sample_ID1": "string",
    "Sample_ID2": "string",
    "PI_HAT": "float",
    "concordance": "float",
    "is_ge_pi_hat": "boolean",
    "is_ge_concordance": "boolean",
}


def read_known_sample_concordance(filename: PathLike) -> pd.DataFrame:
    """Read the known replicate concordance table.

    Returns:
        A table with:

        - ``Subject_ID``
        - ``Sample_ID1``
        - ``Sample_ID2``
        - ``PI_HAT``
        - ``concordance``
        - ``is_ge_pi_hat``
        - ``is_ge_concordance``
    """
    return pd.read_csv(filename, dtype=KNOWN_DTYPES)


def read_unknown_sample_concordance(filename: PathLike) -> pd.DataFrame:
    """Read the unknown replicate concordance table.

    Returns:
        A table with:

        - ``Subject_ID1``
        - ``Subject_ID1``
        - ``Sample_ID1``
        - ``Sample_ID2``
        - ``PI_HAT``
        - ``concordance``
        - ``is_ge_pi_hat``
        - ``is_ge_concordance``
    """
    return pd.read_csv(filename, dtype=UNKNOWN_DTYPES)


@app.command()
def main(
    sample_sheet_csv: Path,
    concordance_csv: Path,
    known_csv: Path,
    known_qc_csv: Path,
    known_study_csv: Path,
    unknown_csv: Path,
):
    # Load sample metadata
    ss = sample_sheet.read(sample_sheet_csv)
    known_replicates = _get_known_replicates(ss)
    sample_to_subject_id = ss.set_index("Sample_ID").Group_By_Subject_ID.rename("Subject_ID")

    # Load sample level concordance information
    concordance = concordance_table.read(concordance_csv).rename(
        {"ID1": "Sample_ID1", "ID2": "Sample_ID2"}, axis=1
    )

    # Get sample concordance for known and unknown replicates
    known_df = _known_replicates_df(concordance, known_replicates, sample_to_subject_id)
    unknown_df = _unknown_replicates_df(concordance, known_replicates, sample_to_subject_id)

    # Save full tables
    known_df.to_csv(known_csv, index=False)
    unknown_df.to_csv(unknown_csv, index=False)

    # Split known replicates into internal controls and study samples
    known_qc, known_study = _split_into_qc_and_study_samples(ss, known_df)
    known_qc.to_csv(known_qc_csv, index=False)
    known_study.to_csv(known_study_csv, index=False)


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


def _known_replicates_df(
    concordance: pd.DataFrame, known_replicates: List, sample_to_subject_id: pd.Series
) -> pd.DataFrame:
    """Filter concordance table for known replicates.

    Pulls known replicate samples from the concordance table. If a replicate
    pair is missing they are filled in as NA.
    """
    if len(known_replicates) == 0:
        # No known replicates
        return pd.DataFrame([], columns=KNOWN_DTYPES.keys())

    res = []
    for s1, s2 in known_replicates:
        query = concordance.query("Sample_ID1 == @s1 & Sample_ID2 == @s2")
        if query.shape[0] > 0:
            res.append(query)
        else:
            # The known replicate pair is missing. This could be due to:
            #   - a sample being removed during call rate filtering
            #   - a pair being pre-filtered because of low concordance (i.e.,
            #     `software_params.ibd_pi_hat_min` [default <= 0.05])
            #   - a pair being pre-filtered because of high concordance (i.e.,
            #     `software_params.ibd_pi_hat_max` [default >= 1.0]). This is less
            #     likely b/c the default is at the upper bound of `PI_HAT`.
            res.append(pd.DataFrame({"Sample_ID1": [s1], "Sample_ID2": [s2]}))

    return (
        pd.concat(res, ignore_index=True)
        .merge(sample_to_subject_id, left_on="Sample_ID1", right_index=True)
        .reset_index()
        .reindex(KNOWN_DTYPES.keys(), axis=1)
    )


def _unknown_replicates_df(
    concordance: pd.DataFrame, known_replicates: List, sample_to_subject_id: pd.Series
) -> pd.DataFrame:
    """Filter concordance table for unknown replicates.

    Pulls unknown replicate samples from the concordance table. These are
    samples that show high levels of concordance but have different
    Subject_IDs.
    """

    def _check(x: pd.Series) -> bool:
        """Test sample pairs.

        Identify sample pairs that are greater than concordance threshold and
        not a known replicate.
        """
        s1, s2 = x.Sample_ID1, x.Sample_ID2
        return x.is_ge_concordance & ((s1, s2) not in known_replicates)

    mask = concordance.apply(_check, axis=1)
    unknown_replicates = concordance[mask]

    return (
        unknown_replicates.merge(
            sample_to_subject_id.rename("Subject_ID1"), left_on="Sample_ID1", right_index=True
        )
        .merge(sample_to_subject_id.rename("Subject_ID2"), left_on="Sample_ID2", right_index=True)
        .reindex(UNKNOWN_DTYPES.keys(), axis=1)
    )


def _split_into_qc_and_study_samples(
    ss: pd.DataFrame, df: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    qc_samples = ss.Sample_ID[ss.is_internal_control].tolist()  # type: ignore # noqa
    return df.query("Sample_ID1 in @qc_samples"), df.query("Sample_ID1 not in @qc_samples")


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {"subject_id_override": None}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({k: Path(v) for k, v in snakemake.output.items()})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
