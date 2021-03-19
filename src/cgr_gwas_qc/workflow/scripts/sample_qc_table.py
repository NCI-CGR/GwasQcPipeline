#!/usr/bin/env python
"""Generate the internal sample QC report.

This script takes various outputs from the GwasQcPipeline and
aggregates/summarizes results into an internal QC report. Here we use a
functional approach to build up the QC report. Each file has it's own
function to parse/summarize its content. This makes it easier to add new
components or change the behavior. Search for `TO-ADD` comments for where you
would have to make modifications to add new components.
"""
import os
from pathlib import Path
from typing import Optional, Sequence
from warnings import warn

import pandas as pd
import typer

from cgr_gwas_qc import load_config
from cgr_gwas_qc.models.config.user_files import Idat
from cgr_gwas_qc.parsers import plink
from cgr_gwas_qc.reporting import CASE_CONTROL_DTYPE, SEX_DTYPE
from cgr_gwas_qc.validators import check_file

app = typer.Typer(add_completion=False)


QC_HEADER = {  # Header for main QC table
    # From Sample Sheet
    "Sample_ID": "string",
    "Group_By_Subject_ID": "string",
    "LIMSSample_ID": "string",
    "LIMS_Individual_ID": "string",
    "SR_Subject_ID": "string",
    "PI_Subject_ID": "string",
    "PI_Study_ID": "string",
    "Project": "string",
    "Project-Sample ID": "string",
    # Generated here
    "num_samples_per_subject": "UInt8",
    "is_internal_control": "boolean",
    "case_control": CASE_CONTROL_DTYPE,
    "is_preflight_exclusion": "boolean",
    "idats_exist": "boolean",
    "expected_sex": SEX_DTYPE,
    "predicted_sex": SEX_DTYPE,
    "X_inbreeding_coefficient": "float",
    "IdatIntensity": "float",
    "AFR": "float",
    "EUR": "float",
    "ASN": "float",
    "Ancestry": "category",
    "Contamination_Rate": "float",
    "Call_Rate_Initial": "float",
    "is_cr1_filtered": "boolean",
    "Call_Rate_1": "float",
    "is_cr2_filtered": "boolean",
    "Call_Rate_2": "float",
    "is_call_rate_filtered": "boolean",
    "is_contaminated": "boolean",
    "is_sex_discordant": "boolean",
    "is_replicate_discordant": "boolean",
    "is_unexpected_replicate": "boolean",
    "Count_of_QC_Issue": "UInt8",
    "identifiler_needed": "boolean",
    "Identifiler_Reason": "string",
    "Subject_Representative": "boolean",
    "Subject_Dropped_From_Study": "boolean",
}


QC_SUMMARY_FLAGS = [  # Set of binary flags used for summarizing sample quality
    "is_call_rate_filtered",
    "is_contaminated",
    "is_sex_discordant",
    "is_replicate_discordant",
    "is_unexpected_replicate",
    # TO-ADD: If you create a new summary binary flag you want to include
    # in the count of QC issues then add the column here.
]

IDENTIFILER_FLAGS = {  # Set of binary flags used to determine if we need to run identifiler
    "is_contaminated": "Contaminated",
    "is_sex_discordant": "Sex Discordant",
    "is_replicate_discordant": "Discordant Replicates",
    "is_unexpected_replicate": "Unexpected Replicate",
    # TO-ADD: If you create a new binary flag do determine if you run
    # identifiler.
}


@app.command()
def main(
    imiss_start: Path = typer.Argument(..., help="Path to plink_start/samples.imiss"),
    imiss_cr1: Path = typer.Argument(..., help="Path to plink_filter_call_rate_1/samples.imiss"),
    imiss_cr2: Path = typer.Argument(..., help="Path to plink_filter_call_rate_2/samples.imiss"),
    sexcheck_cr1: Path = typer.Argument(
        ..., help="Path to plink_filter_call_rate_1/samples.sexcheck"
    ),
    ancestry: Path = typer.Argument(..., help="Path to ancestry/graf_ancestry_calls.txt"),
    known_replicates: Path = typer.Argument(..., help=""),
    unknown_replicates: Path = typer.Argument(..., help=""),
    # TO-ADD: new file to include in QC report. You also need to add this file
    # to the rule in `reporting.smk`
    # Optional inputs
    contam: Optional[Path] = typer.Option(
        None, help="Path to sample_filters/agg_contamination_test.csv"
    ),
    intensity: Optional[Path] = typer.Option(
        None, help="Path to sample_filters/agg_median_idat_intensity.csv"
    ),
    # Outputs
    all_samples: Path = typer.Argument(..., help="Path to all_samples_qc.csv"),
):

    cfg = load_config()
    ss = _wrangle_sample_sheet(cfg.ss, cfg.config.workflow_params.expected_sex_col_name)
    Sample_IDs = ss.index

    ################################################################################
    # Build QC Table
    ################################################################################
    sample_qc = (
        pd.concat(
            [
                ss,
                _check_preflight(cfg.config.Sample_IDs_to_remove, Sample_IDs),
                _check_idats_files(ss, cfg.config.user_files.idat_pattern),
                _read_imiss(imiss_start, Sample_IDs, "Call_Rate_Initial"),
                _read_imiss(imiss_cr1, Sample_IDs, "Call_Rate_1"),
                _read_imiss(imiss_cr2, Sample_IDs, "Call_Rate_2"),
                _read_sexcheck_cr1(sexcheck_cr1, ss.expected_sex),
                _read_ancestry(ancestry, Sample_IDs),
                _read_known_replicates(
                    known_replicates, cfg.config.software_params.dup_concordance_cutoff, Sample_IDs
                ),
                _read_unknown_replicates(unknown_replicates, Sample_IDs),
                _read_contam(contam, cfg.config.software_params.contam_threshold, Sample_IDs),
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
    # Add Call Rate Flags
    # NOTE: Each call rate step will drop samples, so missing samples were
    # filtered in the current or previous step(s). For example, a missing sample
    # `Call_Rate_1` indicates the sample did not pass Call Rate 1 filter or
    # was missing from before. If the sample was missing from before than I am
    # setting it as missing (pd.NA). This will help with data provenance and
    # hopefully make it clearer why a sample was removed.
    cri = sample_qc.Call_Rate_Initial.isna()
    cr1 = sample_qc.Call_Rate_1.isna()
    cr2 = sample_qc.Call_Rate_2.isna()

    sample_qc["is_cr1_filtered"] = cr1
    sample_qc.loc[cri, "is_cr1_filtered"] = pd.NA

    sample_qc["is_cr2_filtered"] = cr2
    sample_qc.loc[cri | cr1, "is_cr2_filtered"] = pd.NA

    # Add Call Rate Summary Flag
    # NOTE: This is `True` if filtered in CR1 or CR2. It is `pd.NA` if missing
    # in the initial data set (cri) and otherwise `False`.
    sample_qc["is_call_rate_filtered"] = cr2
    sample_qc.loc[cri, "is_call_rate_filtered"] = pd.NA

    # Count the number of QC issues
    sample_qc["Count_of_QC_Issue"] = sample_qc[QC_SUMMARY_FLAGS].sum(axis=1).astype(int)

    # Add a flag to run identifiler based if any of these columns are True
    sample_qc["identifiler_needed"] = sample_qc[IDENTIFILER_FLAGS].any(axis=1)
    sample_qc["Identifiler_Reason"] = _identifiler_reason(sample_qc, list(IDENTIFILER_FLAGS))

    # Add flag for which samples to keep as subject
    sample_qc["Subject_Representative"] = _find_study_subject_representative(sample_qc)
    sample_qc["Subject_Dropped_From_Study"] = _find_study_subject_with_no_representative(sample_qc)
    ################################################################################
    # Save Output
    ################################################################################
    _save_qc_table(sample_qc, all_samples)


def read_sample_qc(filename: os.PathLike) -> pd.DataFrame:
    """Read the Sample Level QC Table

    Returns:
        pd.DataFrame:
            Assigning specific data types for optimal parsing.
    """
    return pd.read_csv(filename, dtype=QC_HEADER)


def _wrangle_sample_sheet(sample_sheet: pd.DataFrame, expected_sex_col_name: str) -> pd.DataFrame:
    """Identify expected sex column and count number samples per subject.

    Users can specify which column in the `sample_sheet` holds the expected
    sex information (`config.workflow_params.expected_sex_col_name`). Here we
    rename `expected_sex_col_name` to `expected_sex`.

    We also add the summary column with the number of `Sample_ID`s per
    `Group_By_Subject_ID`.

    Returns:
        pd.DataFrame: The full sample sheet with the following adjustments.
            - Sample_ID (pd.Index)
            - num_samples_per_subject (int): Number of samples per
              `Group_By_Subject_ID`.
            - expected_sex (category): M/F/U based on the expected set column
              set in the config.
            - case_control (category): case/control/qc based on the
              Case/Control_Status in the sample sheet.
    """
    df = sample_sheet.copy()

    df["is_internal_control"] = (df.Sample_Group == "sVALD-001").astype("boolean")

    sex_mapper = {"m": "M", "male": "M", "f": "F", "female": "F"}
    df["expected_sex"] = (
        df[expected_sex_col_name]
        .str.lower()
        .map(lambda sex: sex_mapper.get(sex, "U"))
        .astype(SEX_DTYPE)
    )

    case_control_mapper = {cat.lower(): cat for cat in CASE_CONTROL_DTYPE.categories}
    df["case_control"] = (
        df["Case/Control_Status"]
        .str.lower()
        .map(lambda status: case_control_mapper.get(status, pd.NA))
        .astype(CASE_CONTROL_DTYPE)
    )

    # For internal controls use the `Indentifiler_Sex` column as `expected_sex`
    df.loc[df.is_internal_control, "expected_sex"] = df.loc[
        df.is_internal_control, "Identifiler_Sex"
    ]

    # For internal controls set case_control to qc
    df.loc[df.is_internal_control, "case_control"] = case_control_mapper["qc"]

    # Count the number of samples per subject ID and set Sample_ID as index
    return df.merge(
        df.groupby("Group_By_Subject_ID", dropna=False).size().rename("num_samples_per_subject"),
        on="Group_By_Subject_ID",
    ).set_index("Sample_ID")


def _check_preflight(samples_to_remove: Optional[Sequence[str]], Sample_IDs: pd.Index) -> pd.Series:
    """Checks if any samples were flagged during pre-flight checks.

    Pre-flight checks save samples with missing GTC or IDAT files to the
    config file. This adds a column to the sample_qc table to indicate if a
    sample was excluded due to pre-flight checks.
    """
    if samples_to_remove is None:
        return pd.Series(False, index=Sample_IDs, name="is_preflight_exclusion", dtype="boolean")

    return pd.Series(
        Sample_IDs.isin(samples_to_remove),
        index=Sample_IDs,
        name="is_preflight_exclusion",
        dtype="boolean",
    )


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


def _read_known_replicates(
    file_name: Path, dup_concordance_cutoff: float, Sample_IDs: pd.Index
) -> pd.Series:  # noqa
    """Create a flag of known replicates that show low concordance.

    Given a set of samples that are known to be from the same Subject. Flag
    samples that show low concordance with one or more replicates.

    Returns:
        pd.Series:
            - Sample_ID (pd.Index)
            - is_replicate_discordant (bool): True if replicates show
              a concordance below the supplied threshold. Otherwise False.
   """

    discord_Sample_IDs = (
        pd.read_csv(file_name)  # Subject_ID Sample_ID1 Sample_ID2 Concordance PI_HAT
        .rename({"Concordance": "concordance"}, axis=1)  # legacy uses uppercase
        .query("concordance.notna() & concordance < @dup_concordance_cutoff")
        .loc[:, ("Sample_ID1", "Sample_ID2")]
        .melt()
        .value.unique()
    )  # A set of Sample_IDs that were replicates were not concordant.

    sr = pd.Series(False, index=Sample_IDs, name="is_replicate_discordant")
    sr[sr.index.isin(discord_Sample_IDs)] = True

    return sr


def _read_unknown_replicates(file_name: Path, Sample_IDs: pd.Index) -> pd.Series:
    """Create a flag of for samples that are likely to be unknown replicates.

    Given a set of samples that are not known to be from the same Subject but
    show high concordance.

    Returns:
        pd.Series:
            - Sample_ID (pd.Index)
            - is_replicate_discordant (bool): True if replicates show
              a concordance below the supplied threshold. Otherwise False.
   """

    cord_Sample_IDs = (
        pd.read_csv(file_name)  # Subject_ID1 Subject_ID2 Sample_ID1 Sample_ID2 Concordance PI_HAT
        .loc[:, ("Sample_ID1", "Sample_ID2")]
        .melt()
        .value.unique()
    )  # A set of Sample_IDs that look like a replicate with another sample from a different subject.

    sr = pd.Series(False, index=Sample_IDs, name="is_unexpected_replicate")
    sr[sr.index.isin(cord_Sample_IDs)] = True

    return sr


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


def _check_idats_files(sample_sheet: pd.DataFrame, idat_pattern: Optional[Idat]) -> pd.Series:
    """Check that red and green IDAT files exist.

    Args:
        df: A sample table with at least `Sample_ID` and columns needed to fill wildcards.
        red: A wildcard pattern for red files. Wildcards must be in `df`.
        green: A wildcard pattern for green files. Wildcards must be in `df`.

    Returns:
        pd.Series:
            - Sample_ID (pd.Index)
            - idats_exist (bool): True if both the red and green Idat files existed
    """
    if not idat_pattern:
        # No Idat path specified in config, return all NaN.
        return pd.Series(index=sample_sheet.index, dtype="boolean", name="idats_exist")

    results = []
    for record in sample_sheet.itertuples():
        Sample_ID = record.Index
        red = idat_pattern.red.format(**record._asdict())
        green = idat_pattern.green.format(**record._asdict())
        try:
            check_file(Path(red))
            check_file(Path(green))
            results.append((Sample_ID, True))
        except (FileNotFoundError, PermissionError):
            results.append((Sample_ID, False))

    return pd.Series(dict(results), dtype="boolean", name="idats_exist").rename_axis("Sample_ID")


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


def _find_study_subject_representative(sample_qc: pd.DataFrame) -> pd.Series:
    """Flag indicating which sample to use as subject representative.

    We use a single representative sample for subject level analysis. First
    we remove all internal controls and poor quality samples (is_call_rate_filtered,
    is_contaminated, Replicate Discordance). For subject IDs with multiple
    remaining samples, we select the sample that has the highest Call Rate 2.

    Returns:
        pd.Series:
            - Sample_ID (pd.Index): sorted by `sample_qc.index`
            - A boolean flag where True indicates a sample was used as the
              subject representative.
    """
    return (
        sample_qc.fillna({k: False for k in QC_SUMMARY_FLAGS})  # query breaks if there are NaNs
        .query(
            "not is_internal_control & not is_contaminated & not is_call_rate_filtered & not `is_replicate_discordant`"
        )
        .groupby("Group_By_Subject_ID")  # Group sample by subject id
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
        .Subject_Representative.sum()
        == 0
    ).pipe(lambda x: x.index[x])
    return sample_qc.Group_By_Subject_ID.isin(subject_w_no_rep)


def _save_qc_table(sample_qc: pd.DataFrame, file_name: Path) -> None:
    """Save main QC table."""
    sample_qc.reindex(QC_HEADER, axis=1).to_csv(file_name, index=False)


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {"contam": None, "intensity": None}
        defaults.update({k: (Path(v) if v else None) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({k: Path(v) for k, v in snakemake.output.items()})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
