#!/usr/bin/env python
"""
Split Sample Concordance Tables
-------------------------------

``workflow/scripts/split_sample_concordance.py``

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
    is_internal_control1, boolean, True if the first sample in the comparison is an internal control.
    is_internal_control2, boolean, True if the second sample in the comparison is an internal control.
    is_discordant_replicate, boolean, True if the pair is a known replicate but PLINK, GRAF, and KING could not call it as a replicate
    PLINK_PI_HAT, float, Proportion IBD i.e. ``P(IBD=2) + 0.5*P(IBD=1)``
    PLINK_concordance, float, Proportion IBS2 ``IBS2 / (IBS0 + IBS1 + IBS2)``
    PLINK_is_ge_pi_hat, boolean, True if PI_HAT was greater than ``software_params.pi_hat_cutoff``
    PLINK_is_ge_concordance, boolean, True if concordance was greater than ``software_params.dup_concordance_cutoff``
    GRAF_HGMR, float, Homozygous Genotype Mismatch Rate (%)
    GRAF_AGMR, float, All Genotype Mismatch Rate (%)
    GRAF_relationship, string, relationship determined by sample genotypes.
    KING_Kinship, float, Estimated kinship coefficient from the SNP data
    KING_relationship, string, The assigned relationship based on Kinship

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
    is_internal_control1, boolean, True if the first sample in the comparison is an internal control.
    is_internal_control2, boolean, True if the second sample in the comparison is an internal control.
    PLINK_PI_HAT, float, Proportion IBD i.e. ``P(IBD=2) + 0.5*P(IBD=1)``
    PLINK_concordance, float, Proportion IBS2 ``IBS2 / (IBS0 + IBS1 + IBS2)``
    PLINK_is_ge_pi_hat, boolean, True if PI_HAT was greater than ``software_params.pi_hat_cutoff``
    PLINK_is_ge_concordance, boolean, True if concordance was greater than ``software_params.dup_concordance_cutoff``
    GRAF_HGMR, float, Homozygous Genotype Mismatch Rate (%)
    GRAF_AGMR, float, All Genotype Mismatch Rate (%)
    GRAF_relationship, string, relationship determined by sample genotypes.
    KING_Kinship, float, Estimated kinship coefficient from the SNP data
    KING_relationship, string, The assigned relationship based on Kinship

References:
    - :mod:`cgr_gwas_qc.workflow.scripts.sample_concordance`
"""
from pathlib import Path

import pandas as pd
import typer

from cgr_gwas_qc.typing import PathLike
from cgr_gwas_qc.workflow.scripts import sample_concordance

app = typer.Typer(add_completion=False)

KNOWN_DTYPES = {
    "Subject_ID": "string",
    "Sample_ID1": "string",
    "Sample_ID2": "string",
    "is_internal_control1": "boolean",
    "is_internal_control2": "boolean",
    "is_discordant_replicate": "boolean",
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

UNKNOWN_DTYPES = {
    "Subject_ID1": "string",
    "Subject_ID2": "string",
    "Sample_ID1": "string",
    "Sample_ID2": "string",
    "is_internal_control1": "boolean",
    "is_internal_control2": "boolean",
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


def read_known_sample_concordance(filename: PathLike) -> pd.DataFrame:
    """Read the known replicate concordance table.

    Returns:
        pd.DataFrame
        - Subject_ID
        - Sample_ID1
        - Sample_ID2
        - is_internal_control1
        - is_internal_control2
        - is_discordant_replicate
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
    return pd.read_csv(filename, dtype=KNOWN_DTYPES)


def read_unknown_sample_concordance(filename: PathLike) -> pd.DataFrame:
    """Read the unknown replicate concordance table.

    Returns:
        pd.DataFrame
        - Subject_ID1
        - Subject_ID1
        - Sample_ID1
        - Sample_ID2
        - is_internal_control1
        - is_internal_control2
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
    return pd.read_csv(filename, dtype=UNKNOWN_DTYPES)


@app.command()
def main(
    sample_concordance_csv: Path,
    known_csv: Path,
    known_qc_csv: Path,
    known_study_csv: Path,
    unknown_csv: Path,
):
    # Load sample level concordance information
    concordance = sample_concordance.read(sample_concordance_csv)

    # Save Known Replicates
    known_df = concordance.query("is_expected_replicate").rename(
        {"Subject_ID1": "Subject_ID"}, axis=1
    )
    known_df.reindex(KNOWN_DTYPES, axis=1).to_csv(known_csv, index=False)

    # Known QC only
    known_df.query("is_internal_control1 & is_internal_control2").reindex(
        KNOWN_DTYPES, axis=1
    ).to_csv(known_qc_csv, index=False)

    # Known Study only
    known_df.query("~is_internal_control1 & ~is_internal_control2").reindex(
        KNOWN_DTYPES, axis=1
    ).to_csv(known_study_csv, index=False)

    # Save Unexpected Replicates
    (
        concordance.query("is_unexpected_replicate")
        .reindex(UNKNOWN_DTYPES, axis=1)
        .to_csv(unknown_csv, index=False)
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {
            **{"sample_concordance_csv": Path(snakemake.input[0])},  # type: ignore # noqa
            **{k: Path(v) for k, v in snakemake.output.items()},  # type: ignore # noqa
        }
        main(**defaults)
    else:
        app()
