#!/usr/bin/env python
"""

sample_qc_table.py
------------------

Generate the internal sample QC report.

This script takes various outputs from the GwasQcPipeline and
aggregates/summarizes results into an internal QC report. Here we use a
functional approach to build up the QC report. Each file has it's own
function to parse/summarize its content. This makes it easier to add new
components or change the behavior. Search for `TO-ADD` comments for where you
would have to make modifications to add new components.
"""
from pathlib import Path
from typing import Optional, Sequence
from warnings import warn

import pandas as pd
import typer

from cgr_gwas_qc.parsers import plink, sample_sheet
from cgr_gwas_qc.reporting import CASE_CONTROL_DTYPE, SEX_DTYPE
from cgr_gwas_qc.typing import PathLike
from cgr_gwas_qc.workflow.scripts import sample_concordance
from cgr_gwas_qc.workflow.scripts.snp_qc_table import add_call_rate_flags

app = typer.Typer(add_completion=False)


DTYPES = {  # Header for main QC table
    "Sample_ID": "string",
    "Group_By_Subject_ID": "string",
    "num_samples_per_subject": "UInt8",
    "is_subject_representative": "boolean",
    "subject_dropped_from_study": "boolean",
    "case_control": CASE_CONTROL_DTYPE,
    "is_internal_control": "boolean",
    "is_sample_exclusion": "boolean",
    "idats_exist": "boolean",
    "Call_Rate_Initial": "float",
    "Call_Rate_1": "float",
    "Call_Rate_2": "float",
    "is_cr1_filtered": "boolean",
    "is_cr2_filtered": "boolean",
    "is_call_rate_filtered": "boolean",
    "IdatIntensity": "float",
    "Contamination_Rate": "float",
    "is_contaminated": "boolean",
    "replicate_ids": "string",
    "is_discordant_replicate": "boolean",
    "is_unexpected_replicate": "boolean",
    "expected_sex": SEX_DTYPE,
    "predicted_sex": SEX_DTYPE,
    "X_inbreeding_coefficient": "float",
    "is_sex_discordant": "boolean",
    "AFR": "float",
    "EUR": "float",
    "ASN": "float",
    "Ancestry": "category",
    "Count_of_QC_Issue": "UInt8",
    "is_pass_sample_qc": "boolean",
    "identifiler_needed": "boolean",
    "identifiler_reason": "string",
}


QC_SUMMARY_FLAGS = [  # Set of binary flags used for summarizing sample quality
    "is_call_rate_filtered",
    "is_contaminated",
    "is_sex_discordant",
    "is_discordant_replicate",
    "is_unexpected_replicate",
    # TO-ADD: If you create a new summary binary flag you want to include
    # in the count of QC issues then add the column here.
]

IDENTIFILER_FLAGS = {  # Set of binary flags used to determine if we need to run identifiler
    "is_contaminated": "Contaminated",
    "is_sex_discordant": "Sex Discordant",
    "is_discordant_replicate": "Discordant Replicates",
    "is_unexpected_replicate": "Unexpected Replicate",
    # TO-ADD: If you create a new binary flag do determine if you run
    # identifiler.
}


@app.command()
def main(
    sample_sheet_csv: Path = typer.Argument(..., help="Path to the sample sheet"),
    imiss_start: Path = typer.Argument(..., help="Path to plink_start/samples.imiss"),
    imiss_cr1: Path = typer.Argument(..., help="Path to plink_filter_call_rate_1/samples.imiss"),
    imiss_cr2: Path = typer.Argument(..., help="Path to plink_filter_call_rate_2/samples.imiss"),
    sexcheck_cr1: Path = typer.Argument(
        ..., help="Path to plink_filter_call_rate_1/samples.sexcheck"
    ),
    ancestry: Path = typer.Argument(..., help="Path to ancestry/graf_ancestry_calls.txt"),
    sample_concordance_csv: Path = typer.Argument(..., help=""),
    # Optional inputs
    contam: Optional[Path] = typer.Option(
        None, help="Path to sample_filters/agg_contamination_test.csv"
    ),
    intensity: Optional[Path] = typer.Option(
        None, help="Path to sample_filters/agg_median_idat_intensity.csv"
    ),
    # Params
    contam_threshold: float = typer.Option(..., help="Threshold for contamination."),
    # Outputs
    outfile: Path = typer.Argument(..., help="Path to output csv"),
):

    ss = sample_sheet.read(sample_sheet_csv, remove_exclusions=False).set_index("Sample_ID")
    Sample_IDs = ss.index

    ################################################################################
    # Build QC Table
    ################################################################################
    sample_qc = (
        pd.concat(
            [
                ss,
                _read_imiss(imiss_start, Sample_IDs, "Call_Rate_Initial"),
                _read_imiss(imiss_cr1, Sample_IDs, "Call_Rate_1"),
                _read_imiss(imiss_cr2, Sample_IDs, "Call_Rate_2"),
                _read_sexcheck_cr1(sexcheck_cr1, ss.expected_sex),
                _read_ancestry(ancestry, Sample_IDs),
                _read_concordance(sample_concordance_csv, Sample_IDs),
                _read_contam(contam, contam_threshold, Sample_IDs),
                _read_intensity(intensity, Sample_IDs),
                # TO-ADD: call function you created to parse/summarize new file
            ],
            axis=1,
        )
        .rename_axis("Sample_ID")
        .reset_index()
    )

    ################################################################################
    # Add Summary Columns
    ################################################################################
    add_call_rate_flags(sample_qc)

    # Count the number of QC issues
    sample_qc["Count_of_QC_Issue"] = sample_qc[QC_SUMMARY_FLAGS].sum(axis=1).astype(int)

    # Add a flag to run identifiler based if any of these columns are True
    sample_qc["identifiler_needed"] = sample_qc[IDENTIFILER_FLAGS].any(axis=1)
    sample_qc["identifiler_reason"] = _identifiler_reason(sample_qc, list(IDENTIFILER_FLAGS))

    # Add flag for which samples to keep as subject
    sample_qc["is_pass_sample_qc"] = _check_pass_qc(sample_qc)
    sample_qc["is_subject_representative"] = _find_study_subject_representative(sample_qc)
    sample_qc["subject_dropped_from_study"] = _find_study_subject_with_no_representative(sample_qc)
    ################################################################################
    # Save Output
    ################################################################################
    _save_qc_table(sample_qc, outfile)


def read(filename: PathLike) -> pd.DataFrame:
    """Read the Sample Level QC Table

    Returns:
        pd.DataFrame:
            Assigning specific data types for optimal parsing.
    """
    return pd.read_csv(filename, dtype=DTYPES)


def _read_imiss(filename: Path, Sample_IDs: pd.Index, col_name: str) -> pd.Series:
    """Read the starting call rates.

    Returns:
        pd.Series:
            - Sample_ID (pd.Index)
            - `col_name` (float): The starting call rate calculated as `1 - F_MISS`.
    """
    return (
        plink.read_imiss(filename)
        .rename_axis("Sample_ID")
        .assign(**{col_name: lambda x: 1 - x.F_MISS})
        .reindex(Sample_IDs)
        .loc[:, col_name]
    )


def _read_sexcheck_cr1(filename: Path, expected_sex: pd.Series) -> pd.DataFrame:
    """Read sex predictions and summarize.

    Read PLINK sex prediction file. Convert the `predicted_sex` indicator
    variable to M/F designations. Compare predicted results with the expected
    sexes and create a summary column `is_sex_discordant` if predicted/expected sex
    calls do not match.

    Returns:
        pd.DataFrame:
            - Sample_ID (pd.Index)
            - X_inbreeding_coefficient (float64): PLINK's inbreeding coefficient
              from sexcheck.
            - predicted_sex (str): M/F/U based on PLINK sex predictions.
              are different. U if prediction was U.
            - is_sex_discordant (bool): True if SexMatch == "N"
    """
    plink_sex_code = {0: "U", 1: "M", 2: "F"}
    df = (
        plink.read_sexcheck(filename)
        .rename_axis("Sample_ID")
        .rename({"F": "X_inbreeding_coefficient"}, axis=1)
        .assign(predicted_sex=lambda x: x.SNPSEX.map(plink_sex_code))
        .astype({"predicted_sex": SEX_DTYPE})
        .reindex(expected_sex.index)
        .reindex(["X_inbreeding_coefficient", "predicted_sex"], axis=1)
    )

    # Update PLINK predicted_sex Calls
    # TODO: Decide if we want to keep this logic from the legacy workflow. See
    # http://10.133.130.114/jfear/GwasQcPipeline/issues/35
    df.loc[df.X_inbreeding_coefficient < 0.5, "predicted_sex"] = "F"
    df.loc[df.X_inbreeding_coefficient >= 0.5, "predicted_sex"] = "M"
    df["is_sex_discordant"] = (df.predicted_sex != expected_sex).astype("boolean")
    df.loc[
        df.X_inbreeding_coefficient.isnull() | (df.predicted_sex == "U"), "is_sex_discordant"
    ] = pd.NA  # If we could not predict sex then label as U

    return df


def _read_ancestry(file_name: Path, Sample_IDs: pd.Index) -> pd.DataFrame:
    """Read ancestry calls from GRAF or SNPweights.

    Raises:
        DeprecationWarning: SNPweights parsing is only maintained for testing
          purposes please use GRAF.

    """
    if "P_f (%)" not in file_name.read_text():
        warn("Please use GRAF instead of SNPweights.", DeprecationWarning)
        return _read_SNPweights(file_name, Sample_IDs)

    return _read_GRAF(file_name, Sample_IDs)


def _read_GRAF(file_name: Path, Sample_IDs: pd.Index) -> pd.DataFrame:
    """Read GRAF's ancestry calls.

    Returns:
        pd.DataFrame:
            - Sample_ID (pd.Index)
            - AFR (float): The proportion African ancestry (`P_f(%) / 100`).
            - EUR (float): The proportion African ancestry (`P_e(%) / 100`).
            - ASN (float): The proportion African ancestry (`P_a(%) / 100`).
            - Ancestry (str): The ancestry label assigned by GRAF. See tables
              2 and 3 in their manuscript_ for the exact criteria.

    .. _manuscript: https://pubmed.ncbi.nlm.nih.gov/31151998/

    """
    return (
        pd.read_csv(file_name, sep="\t")
        .rename({"Subject": "Sample_ID"}, axis=1)
        .assign(Ancestry=lambda x: x["Computed population"].str.replace(" ", "_"))
        .assign(AFR=lambda x: x["P_f (%)"] / 100)
        .assign(EUR=lambda x: x["P_e (%)"] / 100)
        .assign(ASN=lambda x: x["P_a (%)"] / 100)
        .set_index("Sample_ID")
        .loc[:, ("AFR", "EUR", "ASN", "Ancestry")]
        .reindex(Sample_IDs)
    )


def _read_SNPweights(file_name: Path, Sample_IDs: pd.Index) -> pd.DataFrame:
    """Read SNPweight's ancestry calls.

    Returns:
        pd.DataFrame:
            - Sample_ID (pd.Index)
            - AFR (float): The proportion African ancestry.
            - EUR (float): The proportion African ancestry.
            - ASN (float): The proportion African ancestry.
            - Ancestry (str): The assigned ancestry label.

    .. _manuscript: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3661048/

    """
    return (
        pd.read_csv(file_name)
        .rename({"ID": "Sample_ID"}, axis=1)
        .set_index("Sample_ID")
        .loc[:, ("AFR", "EUR", "ASN", "Ancestry")]
        .reindex(Sample_IDs)
    )


def _read_concordance(filename: Path, Sample_IDs: pd.Index) -> pd.Series:  # noqa
    """Create a flag of known replicates that show low concordance.

    Given a set of samples that are known to be from the same Subject. Flag
    samples that show low concordance with one or more replicates.

    Returns:
        pd.Series:
            - Sample_ID (pd.Index)
            - is_discordant_replicate (bool): True if replicates show
              a concordance below the supplied threshold. Otherwise False.
            - is_unexpected_replicate (bool): True if replicates are from two
              different subjects.
   """
    return (
        sample_concordance.read(filename)
        .melt(
            id_vars=["is_discordant_replicate", "is_unexpected_replicate"],
            value_vars=["Sample_ID1", "Sample_ID2"],
            var_name="To_Drop",
            value_name="Sample_ID",
        )
        .drop("To_Drop", axis=1)
        .groupby("Sample_ID")
        .max()  # Flag a sample as True if it is True for any comparison.
        .reindex(Sample_IDs)
    )


def _read_contam(
    file_name: Optional[Path], contam_threshold: float, Sample_IDs: pd.Index
) -> pd.DataFrame:
    """Parse verifyIDintensity contamination information.

    Returns:
        pd.DataFrame:
            - Sample_ID (pd.index)
            - Contamination_Rate (float): The contamination rate as estimated
              by verifyIDintensity.
            - is_contaminated (bool): True if the contamination rate is greater
              than the supplied threshold.
    """

    if file_name is None:
        return pd.DataFrame(index=Sample_IDs, columns=["Contamination_Rate", "is_contaminated"])

    df = (
        pd.read_csv(file_name)
        .rename({"ID": "Sample_ID", "%Mix": "Contamination_Rate"}, axis=1)
        .set_index("Sample_ID")
    )

    df["is_contaminated"] = df.Contamination_Rate > contam_threshold
    df.loc[df.Contamination_Rate.isna(), "is_contaminated"] = False

    return df.reindex(Sample_IDs)[["Contamination_Rate", "is_contaminated"]]


def _read_intensity(file_name: Optional[Path], Sample_IDs: pd.Index) -> pd.Series:
    """Parse the median Idat intensity table.

    Returns:
        pd.Series:
            - Sample_ID (pd.Index)
            - IdatIntensity (float): The median Idat intensity.
    """
    if file_name is None:
        return pd.Series(index=Sample_IDs, dtype="float", name="IdatIntensity")

    return (
        pd.read_csv(file_name)
        .rename(
            {
                "SampId": "Sample_ID",
                "median_intensity": "IdatIntensity",
                "MedianIntensity": "IdatIntensity",  # legacy uses different name
            },
            axis=1,
        )
        .set_index("Sample_ID")
        .reindex(Sample_IDs)
        .IdatIntensity
    )


# TO-ADD: Add a parsing/summary function that returns a Series or DataFrame indexed by Sample_ID


def _identifiler_reason(sample_qc: pd.DataFrame, cols: Sequence[str]):
    """Summary string of the reason for needing identifiler.

    If `identifiler_needed` then, if the binary flag in `cols` is True, then
    concatenate the column names.

    Example:
        >>> cols = ["is_sex_discordant", "is_contaminated"]
        >>> df.values == np.ndarray([[True, True], [True, False], [False, False]])
        >>> _identifiler_reason(df, cols)
        pd.Series(["is_sex_discordant;is_contaminated", "is_sex_discordant", ""])
    """

    def reason_string(row: pd.Series) -> str:
        if row.identifiler_needed:
            flags = row[cols].fillna(False)
            return ";".join(IDENTIFILER_FLAGS.get(x, x) for x in flags.index[flags])
        return ""

    return sample_qc.apply(reason_string, axis=1)


def _check_pass_qc(sample_qc: pd.DataFrame) -> pd.Series:
    """True if a sample passed on sample level QC checks.

    We remove all internal controls (``is_internal_control``) and poor
    quality samples (``is_call_rate_filtered``, ``is_contaminated``,
    ``is_discordant_replicate``).
    """
    df = sample_qc.copy()
    df.fillna({k: False for k in QC_SUMMARY_FLAGS}, inplace=True)  # query breaks if there are NaNs
    return (
        ~df.is_internal_control
        & ~df.is_contaminated
        & ~df.is_call_rate_filtered
        & ~df.is_discordant_replicate
    )


def _find_study_subject_representative(sample_qc: pd.DataFrame) -> pd.Series:
    """Flag indicating which sample to use as subject representative.

    We use a single representative sample for subject level analysis. First
    we keep all samples (``is_pass_sample_qc``). For subject IDs with
    multiple remaining samples, we select the sample that has the highest
    Call Rate 2.

    Returns:
        pd.Series:
            - Sample_ID (pd.Index): sorted by `sample_qc.index`
            - A boolean flag where True indicates a sample was used as the
              subject representative.
    """
    return (
        sample_qc.query("is_pass_sample_qc")
        .groupby("Group_By_Subject_ID")
        .apply(
            lambda x: x.Call_Rate_2 == x.Call_Rate_2.max()
        )  # Select the sample with highest call rate as representative
        .droplevel(0)  # drop the subject label b/c don't need
        .reindex(
            sample_qc.index
        )  # Add the samples that were filtered by the query step and make sure everything aligns
        .fillna(False)
    )


def _find_study_subject_with_no_representative(sample_qc: pd.DataFrame) -> pd.Series:
    """Flag indicating which subjects have representative sample.

    This flag excludes internal controls which by nature are ignored in
    subject level analysis.

    Returns:
        pd.Series:
            - Sample_ID (pd.Index): sorted by `sample_qc.index`
            - A boolean flag where True indicates the subject (i.e., all
              samples) has no representative sample.
    """
    subject_w_no_rep = (
        sample_qc.query("not is_internal_control")
        .groupby("Group_By_Subject_ID")
        .is_subject_representative.sum()
        == 0
    ).pipe(lambda x: x.index[x])
    return sample_qc.Group_By_Subject_ID.isin(subject_w_no_rep)


def _save_qc_table(sample_qc: pd.DataFrame, file_name: Path) -> None:
    """Save main QC table."""
    sample_qc.reindex(DTYPES, axis=1).to_csv(file_name, index=False)


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {"contam": None, "intensity": None}
        defaults.update({k: (Path(v) if v else None) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({k: v for k, v in snakemake.params.items()})  # type: ignore # noqa
        defaults.update({"outfile": Path(snakemake.output[0])})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
