from pathlib import Path
from typing import Optional

import pandas as pd
import typer

from cgr_gwas_qc.parsers.sample_sheet import SampleSheet

app = typer.Typer(add_completion=False)

minConcordance = 1.0


@app.command()
def main(
    sample_sheet: Path = typer.Argument(..., help="Path to the sample sheet CSV."),
    imiss_file: Path = typer.Argument(
        ..., help="Path to the `plink_filter_call_rate_2/samples.imiss`"
    ),
    ibd_file: Path = typer.Argument(..., help="Path to the `ibd/samples.genome` file."),
    dup_concordance_cutoff: float = typer.Argument(
        ..., help="The concordance threshold for samples from different subjects."
    ),
    known_file: Path = typer.Argument(..., help="Path to save known concordant samples."),
    unknown_file: Path = typer.Argument(..., help="Path to save unknown condorant samples."),
    subject_id_override: Optional[str] = typer.Argument(
        None,
        help="Specify which column corresponds to subjects. [Deprecated: Use the `Group_By` column sample sheet]",
    ),
):
    sample_metadata = read_sample_metadata(sample_sheet, imiss_file, subject_id_override)
    sample_concordance = read_pairwise_sample_concordance(ibd_file)

    # Merge on sample metadata for Sample 1 and Sample 2 in the concordance
    # comparison. Add suffix to metadata columns.
    df = sample_concordance.merge(
        sample_metadata.add_suffix("1"), on="Sample_ID1", how="left"
    ).merge(sample_metadata.add_suffix("2"), on="Sample_ID2", how="left")

    # Save concordance information for samples from the same subject
    create_known_concordant_table(df).to_csv(known_file, index=False)

    # Save a list of samples from different subjects but with high concordance
    create_unknown_concordant_table(df, dup_concordance_cutoff).to_csv(unknown_file, index=False)


def read_sample_metadata(
    sample_sheet: Path, call_rate: Path, subject_id_override: Optional[str] = None
) -> pd.DataFrame:
    """Read in sample metadata from sample sheet and call rates.

    Returns:
        DataFrame with ["Sample_ID", "Subject_ID", "call_rate"]
    """
    return (
        SampleSheet(sample_sheet)
        .add_group_by_column(subject_id_override)
        .data.merge(read_imiss_file(call_rate), how="left")
        .rename({"Group_By_Subject_ID": "Subject_ID"}, axis=1)
        .loc[:, ("Sample_ID", "Subject_ID", "call_rate")]
    )


def read_pairwise_sample_concordance(file_name: Path) -> pd.DataFrame:
    """Parse IBD file and calculate sample concordance.

    Returns:
        DataFrame with ["Sample_ID1", "Sample_ID2", "PI_HAT", "concordance"].
          Sample concordance is `IBS2 / (IBS0 + IBS1 + IBS2)`.
    """
    return (
        pd.read_csv(file_name, delim_whitespace=True)
        .assign(concordance=lambda x: x.IBS2 / (x.IBS0 + x.IBS1 + x.IBS2))
        .rename({"IID1": "Sample_ID1", "IID2": "Sample_ID2"}, axis=1)
        .loc[:, ("Sample_ID1", "Sample_ID2", "PI_HAT", "concordance")]
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
    """
    return df.query("Subject_ID1 != Subject_ID2 & concordance > @dup_concordance_cutoff").loc[
        :, ("Subject_ID1", "Subject_ID2", "Sample_ID1", "Sample_ID2", "concordance", "PI_HAT")
    ]


if __name__ == "__main__":
    if "snakemake" in locals():
        main(
            snakemake.input.sampleSheet,  # type: ignore # noqa
            snakemake.input.imiss3,  # type: ignore # noqa
            snakemake.input.ibd,  # type: ignore # noqa
            snakemake.params.concordance_threshold,  # type: ignore # noqa
            snakemake.output.known,  # type: ignore # noqa
            snakemake.output.unknown,  # type: ignore # noqa
            snakemake.params.subject_id_override,  # type: ignore # noqa
        )
    else:
        app()
