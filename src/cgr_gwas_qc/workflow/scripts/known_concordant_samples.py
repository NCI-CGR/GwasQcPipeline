#!/usr/bin/env python
"""Create tables of sample concordance.

This script uses sample concordance based on pairwise IBD estimations. We
expect samples from the same subject to show high concordance while samples
from different subjects should not be related. This script outputs 4 tables:

- Known concordant samples (Full Table)
- Known concordant samples (QC Samples Only)
- Known concordant samples (Study Samples Only)
- Unknown concordant samples (Full Table)
"""
from pathlib import Path
from typing import Optional

import pandas as pd
import typer

from cgr_gwas_qc.parsers.sample_sheet import SampleSheet

app = typer.Typer(add_completion=False)


@app.command()
def main(
    sample_sheet: Path = typer.Argument(..., help="Path to the sample sheet CSV."),
    imiss: Path = typer.Argument(..., help="Path to the `plink_filter_call_rate_2/samples.imiss`"),
    concordance: Path = typer.Argument(..., help="Path to the `ibd/samples.genome` file."),
    known: Path = typer.Argument(..., help="Path to save known concordant samples."),
    known_qc: Path = typer.Argument(..., help="Path to save known concordant internal QC samples."),
    known_study: Path = typer.Argument(..., help="Path to save known concordant study samples."),
    unknown: Path = typer.Argument(..., help="Path to save unknown condorant samples."),
    concordance_threshold: float = typer.Argument(
        ..., help="The concordance threshold for samples from different subjects."
    ),
    subject_id_override: Optional[str] = typer.Argument(
        None,
        help="Specify which column corresponds to subjects. [Deprecated: Use the `Group_By` column sample sheet]",
    ),
):
    ################################################################################
    # Pairwise sample concordance + sample metadata
    ################################################################################
    sample_concordance = pd.read_csv(concordance).rename(
        {"IID1": "Sample_ID1", "IID2": "Sample_ID2"}, axis=1
    )
    sample_metadata = read_sample_metadata(sample_sheet, imiss, subject_id_override)

    # Add metadata for Sample 1 and Sample 2 in pairwise table
    df = sample_concordance.merge(
        sample_metadata.add_suffix("1"), on="Sample_ID1", how="left"
    ).merge(sample_metadata.add_suffix("2"), on="Sample_ID2", how="left")

    ################################################################################
    # Save Known Concordant Samples
    # (i.e., samples from the same subject)
    ################################################################################
    # Full Table
    known_df = create_known_concordant_table(df)
    known_df.to_csv(known, index=False)  # Full table

    # QC samples only
    qc_subjects = sample_metadata.query("Sample_Group == 'sVALD-001'").Subject_ID.values  # noqa
    known_df.query("Subject_ID in @qc_subjects").to_csv(known_qc, index=False)

    # Study samples only
    known_df.query("Subject_ID not in @qc_subjects").to_csv(known_study, index=False)

    ################################################################################
    # Save Unknown Concordant Samples
    # (i.e., highly concordant samples from different subjects)
    ################################################################################
    unknown_df = create_unknown_concordant_table(df, concordance_threshold)
    unknown_df.to_csv(unknown, index=False)


def read_sample_metadata(
    sample_sheet: Path, call_rate: Path, subject_id_override: Optional[str] = None
) -> pd.DataFrame:
    """Read in sample metadata from sample sheet and call rates.

    Returns:
        DataFrame with ["Sample_ID", "Subject_ID", "Sample_Group", "call_rate"]
    """
    return (
        SampleSheet(sample_sheet)
        .add_group_by_column(subject_id_override)
        .data.merge(read_imiss_file(call_rate), how="left")
        .rename({"Group_By_Subject_ID": "Subject_ID"}, axis=1)
        .loc[:, ("Sample_ID", "Subject_ID", "Sample_Group", "call_rate")]
    )


def read_imiss_file(file_name: Path) -> pd.Series:
    """Read imiss file and calculate call rate.

    Returns:
        Series (Sample_ID, call_rate) where call_rate is 1 - F_MISS.
    """
    return (
        pd.read_csv(file_name, delim_whitespace=True)
        .assign(call_rate=lambda df: 1.0 - df.F_MISS)
        .rename({"IID": "Sample_ID"}, axis=1)
        .loc[:, ("Sample_ID", "call_rate")]
        .squeeze()
    )


def create_known_concordant_table(df: pd.DataFrame) -> pd.DataFrame:
    """Save samples known to be from the same subject.

    - If missing PI_HAT estimates set defaults (PI_HAT = 0.05, concordance =
      min concordance)
    - If either sample missing call rate then set (PI_HAT and concordance to
      Nan)

    Returns:
        DataFrame with ["Subject_ID", "Sample_ID1", "Sample_ID2",
          "concordance", "PI_HAT"]
    """
    _df = df.copy()

    missing_pi_hat = _df.PI_HAT.isna()
    _df.loc[missing_pi_hat, "PI_HAT"] = 0.05
    _df.loc[missing_pi_hat, "concordance"] = min(1.0, _df.concordance.min())

    missing_call_rate = _df.call_rate1.isna() | _df.call_rate2.isna()
    _df.loc[missing_call_rate, "PI_HAT"] = "NA"
    _df.loc[missing_call_rate, "concordance"] = "NA"

    return (
        _df.query("Subject_ID1 == Subject_ID2")
        .rename({"Subject_ID1": "Subject_ID"}, axis=1)
        .loc[:, ("Subject_ID", "Sample_ID1", "Sample_ID2", "concordance", "PI_HAT")]
    )


def create_unknown_concordant_table(
    df: pd.DataFrame, dup_concordance_cutoff: float
) -> pd.DataFrame:  # noqa
    """Identify concordant samples from different subjects.

    1. Samples have different subject IDs
    2. Sample concordance is > `dup_concordance_cutoff`

    Returns:
        DataFrame with ["Subject_ID1", "Subject_ID2", "Sample_ID1",
        "Sample_ID2", "concordance", "PI_HAT"]
    """
    return df.query("Subject_ID1 != Subject_ID2 & concordance > @dup_concordance_cutoff").loc[
        :, ("Subject_ID1", "Subject_ID2", "Sample_ID1", "Sample_ID2", "concordance", "PI_HAT")
    ]


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {"subject_id_override": None}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({k: Path(v) for k, v in snakemake.output.items()})  # type: ignore # noqa
        defaults.update(snakemake.params)  # type: ignore # noqa
        main(**defaults)
    else:
        app()
