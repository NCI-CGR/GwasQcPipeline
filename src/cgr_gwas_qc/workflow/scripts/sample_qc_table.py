#!/usr/bin/env python
"""
Internal Sample QC Report
-------------------------

**Script**: ``workflow/scripts/sample_qc_table.py``

This script takes various outputs from the GwasQcPipeline and
aggregates/summarizes results into an internal QC report. Here we use a
functional approach to build up the QC report. Each file has it's own
function to parse/summarize its content. This makes it easier to add new
components or change the behavior. Search for `TO-ADD` comments for where you
would have to make modifications to add new components.

**Output**: ``sample_level/sample_qc.csv``

.. csv-table::
    :header: name, dtype, description

    Sample_ID, string, The sample identifier used by the workflow.
    Group_By_Subject_ID, string, The subject identifier used by the workflow
    num_samples_per_subject, UInt32, The total number of ``Sample_ID`` per ``Group_By_Subject_ID``.
    analytic_exclusion, boolean, True if the sample failed QC criteria.
    num_analytic_exclusion, boolean, The number of QC criteria the sample failed.
    analytic_exclusion_reason, string, A list of QC criteria the sample failed.
    is_subject_representative, boolean, True if the ``Sample_ID`` is used to represent the subject.
    subject_dropped_from_study, boolean, True if there are no samples passing QC for this subject.
    case_control, CASE_CONTROL_DTYPE, Phenotype status [``Case`` | ``Control`` | ``QC`` | ``Unknown``].
    is_internal_control, boolean, True if the sample is an internal control.
    is_sample_exclusion, boolean, True if the sample is excluded in the config or missing IDAT/GTC files.
    is_user_exclusion, boolean, True if the sample is excluded in the config
    is_missing_idats, boolean, True if the IDAT files were missing during pre-flight.
    is_missing_gtc, boolean, True if the IDAT files were missing during pre-flight.
    Call_Rate_Initial, float, The Initial sample call rate before any filters.
    Call_Rate_1, float, The sample call rate after the first filter.
    Call_Rate_2, float, The sample call rate after the second filter.
    is_cr1_filtered, boolean, True if a sample was removed by the first call rate filter.
    is_cr2_filtered, boolean, True if a sample was removed by the second call rate filter.
    is_call_rate_filtered, boolean, True if a sample was removed by either call rate filter.
    IdatIntensity, float, The median IDAT intensity.
    Contamination_Rate, float, The SNPweights contamination rate (%Mix)
    is_contaminated, boolean, True if a sample exceed the contamination rate threshold.
    replicate_ids, string, Concatenated Sample_IDs that are replicates.
    is_discordant_replicate, boolean, True if the replicates had low concordance.
    expected_sex, SEX_DTYPE, The expected sex from the provided sample sheet.
    predicted_sex, SEX_DTYPE, The predicted sex based on X chromosome heterozygosity.
    X_inbreeding_coefficient, float, The X chromosome F coefficient.
    is_sex_discordant, boolean, True if expected and predicted sex did not match.
    AFR, float, Proportion African ancestry.
    EUR, float, Proportion European ancestry.
    ASN, float, Proportion East Asian ancestry.
    Ancestry, category, Assigned ancestral population (GRAF).
    identifiler_needed, boolean, True if the lab needs to run Identifiler.
    identifiler_reason, string, The QC problem that needs checked with Identifiler.

"""
from pathlib import Path
from typing import Mapping, Optional
from warnings import warn

import pandas as pd
import typer

from cgr_gwas_qc.parsers import plink, sample_sheet
from cgr_gwas_qc.reporting import CASE_CONTROL_DTYPE, SEX_DTYPE
from cgr_gwas_qc.typing import PathLike
from cgr_gwas_qc.workflow.scripts import agg_contamination, sample_concordance
from cgr_gwas_qc.workflow.scripts.snp_qc_table import add_call_rate_flags

app = typer.Typer(add_completion=False)


DTYPES = {  # Header for main QC table
    "Sample_ID": "string",
    "Group_By_Subject_ID": "string",
    "num_samples_per_subject": "UInt32",
    "analytic_exclusion": "boolean",
    "num_analytic_exclusion": "UInt8",
    "analytic_exclusion_reason": "string",
    "is_subject_representative": "boolean",
    "subject_dropped_from_study": "boolean",
    "case_control": CASE_CONTROL_DTYPE,
    "is_internal_control": "boolean",
    "is_sample_exclusion": "boolean",
    "is_user_exclusion": "boolean",
    "is_missing_idats": "boolean",
    "is_missing_gtc": "boolean",
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
    "expected_sex": SEX_DTYPE,
    "predicted_sex": SEX_DTYPE,
    "X_inbreeding_coefficient": "float",
    "is_sex_discordant": "boolean",
    "AFR": "float",
    "EUR": "float",
    "ASN": "float",
    "Ancestry": "category",
    "identifiler_needed": "boolean",
    "identifiler_reason": "string",
}


def read(filename: PathLike) -> pd.DataFrame:
    """Read the Sample Level QC Table

    Returns:
        pd.DataFrame:
        - Sample_ID
        - Group_By_Subject_ID
        - num_samples_per_subject
        - analytic_exclusion
        - num_analytic_exclusion
        - analytic_exclusion_reason
        - is_subject_representative
        - subject_dropped_from_study
        - case_control
        - is_internal_control
        - is_sample_exclusion
        - is_user_exclusion
        - is_missing_idats
        - is_missing_gtc
        - Call_Rate_Initial
        - Call_Rate_1
        - Call_Rate_2
        - is_cr1_filtered
        - is_cr2_filtered
        - is_call_rate_filtered
        - IdatIntensity
        - Contamination_Rate
        - is_contaminated
        - replicate_ids
        - is_discordant_replicate
        - expected_sex
        - predicted_sex
        - X_inbreeding_coefficient
        - is_sex_discordant
        - AFR
        - EUR
        - ASN
        - Ancestry
        - identifiler_needed
        - identifiler_reason
    """
    return pd.read_csv(filename, dtype=DTYPES)


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
    remove_contam: bool = True,
    remove_rep_discordant: bool = True,
    # Outputs
    outfile: Path = typer.Argument(..., help="Path to output csv"),
):
    ss = sample_sheet.read(sample_sheet_csv, remove_exclusions=False).set_index("Sample_ID")
    sample_qc = build(
        ss,
        imiss_start,
        imiss_cr1,
        imiss_cr2,
        sexcheck_cr1,
        ancestry,
        sample_concordance_csv,
        contam,
        intensity,
    )
    add_qc_columns(
        sample_qc, remove_contam, remove_rep_discordant,
    )
    save(sample_qc, outfile)


def build(
    ss: pd.DataFrame,
    imiss_start: Path,
    imiss_cr1: Path,
    imiss_cr2: Path,
    sexcheck_cr1: Path,
    ancestry: Path,
    sample_concordance_csv: Path,
    contam: Optional[Path],
    intensity: Optional[Path],
) -> pd.DataFrame:
    Sample_IDs = ss.index
    return (
        pd.concat(
            [
                ss,
                _read_imiss(imiss_start, Sample_IDs, "Call_Rate_Initial"),
                _read_imiss(imiss_cr1, Sample_IDs, "Call_Rate_1"),
                _read_imiss(imiss_cr2, Sample_IDs, "Call_Rate_2"),
                _read_sexcheck_cr1(sexcheck_cr1, ss.expected_sex),
                _read_ancestry(ancestry, Sample_IDs),
                _read_concordance(sample_concordance_csv, Sample_IDs),
                _read_contam(contam, Sample_IDs),
                _read_intensity(intensity, Sample_IDs),
                # TO-ADD: call function you created to parse/summarize new file
            ],
            axis=1,
        )
        .rename_axis("Sample_ID")
        .reset_index()
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
        .reindex([col_name], axis=1)
        .squeeze()
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
    df.loc[df.X_inbreeding_coefficient.isnull(), "predicted_sex"] = "U"
    df.loc[df.X_inbreeding_coefficient < 0.5, "predicted_sex"] = "F"
    df.loc[df.X_inbreeding_coefficient >= 0.5, "predicted_sex"] = "M"
    df["is_sex_discordant"] = (df.predicted_sex != expected_sex).astype("boolean")

    # If expected or predicted sex are unknown then set discordance to NA
    df.loc[(expected_sex == "U") | (df.predicted_sex == "U"), "is_sex_discordant"] = pd.NA

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


def _read_concordance(filename: Path, Sample_IDs: pd.Index) -> pd.DataFrame:
    """Create a flag of known replicates that show low concordance.

    Given a set of samples that are known to be from the same Subject. Flag
    samples that show low concordance with one or more replicates.

    Returns:
        pd.Series:
            - Sample_ID (pd.Index)
            - is_discordant_replicate (bool): True if replicates show
              a concordance below the supplied threshold. Otherwise False.
    """
    df = sample_concordance.read(filename)
    return (
        df.melt(
            id_vars=["is_discordant_replicate"],
            value_vars=["Sample_ID1", "Sample_ID2"],
            var_name="To_Drop",
            value_name="Sample_ID",
        )
        .drop("To_Drop", axis=1)
        .groupby("Sample_ID")
        .max()  # Flag a sample as True if it is True for any comparison.
        .astype("boolean")
        .reindex(Sample_IDs)
    )


def _read_contam(file_name: Optional[Path], Sample_IDs: pd.Index) -> pd.DataFrame:
    """Parse verifyIDintensity contamination information.

    Returns:
        pd.DataFrame:
        - Sample_ID (pd.index)
        - Contamination_Rate (float): The contamination rate as estimated
            by verifyIDintensity (%Mix).
        - is_contaminated (bool): True if the contamination rate is greater
            than the supplied threshold.
    """

    if file_name is None:
        return pd.DataFrame(
            index=Sample_IDs, columns=["Contamination_Rate", "is_contaminated"],
        ).astype({"Contamination_Rate": "float", "is_contaminated": "boolean"})

    return (
        agg_contamination.read(file_name)
        .rename({"%Mix": "Contamination_Rate"}, axis=1)
        .set_index("Sample_ID")
        .reindex(["Contamination_Rate", "is_contaminated"], axis=1)
        .reindex(Sample_IDs)
    )


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


def add_qc_columns(
    sample_qc: pd.DataFrame, remove_contam: bool, remove_rep_discordant: bool,
) -> pd.DataFrame:
    add_call_rate_flags(sample_qc)
    _add_identifiler(sample_qc)
    _add_analytic_exclusion(
        sample_qc, remove_contam, remove_rep_discordant,
    )
    _add_subject_representative(sample_qc)
    _add_subject_dropped_from_study(sample_qc)

    return sample_qc


def _add_identifiler(sample_qc: pd.DataFrame) -> pd.DataFrame:
    """Add a flag to run identifiler based if any of these columns are True"""
    identifiler_flags = {  # Set of binary flags used to determine if we need to run identifiler
        "is_contaminated": "Contamination",
        "is_sex_discordant": "Sex Discordance",
        "is_discordant_replicate": "Replicate Discordance",
        # TO-ADD: If you create a new binary flag do determine if you run
        # identifiler.
    }
    sample_qc["identifiler_needed"] = sample_qc[identifiler_flags.keys()].any(axis=1)
    sample_qc["identifiler_reason"] = _get_reason(sample_qc, identifiler_flags)
    return sample_qc


def _get_reason(sample_qc: pd.DataFrame, flags: Mapping[str, str]):
    """Summary string of the reason.

    Given a set of boolean columns, this will create a series with the column
    names where there were True values.

    Example:
        >>> df.columns == ["is_sex_discordant", "is_contaminated"]
        >>> df.values == np.ndarray([[True, True], [True, False], [False, False]])
        >>> _identifiler_reason(df, {"is_sex_discordant": "Sex Discordant", "is_contaminated": "Contaminated"})
        pd.Series(["Sex Discordant;Contaminated", "Sex Discordant", ""])
    """
    cols = list(flags.keys())

    def reason_string(row: pd.Series) -> str:
        row_flags = row.reindex(cols).fillna(False)
        if row_flags.any():
            return "|".join(sorted(flags.get(x, x) for x in row_flags.index[row_flags]))
        return pd.NA

    return sample_qc.apply(reason_string, axis=1)


def _add_analytic_exclusion(
    sample_qc: pd.DataFrame, remove_contam: bool, remove_rep_discordant: bool,
) -> pd.DataFrame:
    """Adds a flag to remove samples based on provided conditions.

    Some of the conditions are user configurable, but we always exclude
    samples with low call rate and samples marked for exclusion during
    pre-flight checks.
    """
    exclusion_criteria = {
        "is_user_exclusion": "User Exclusion",
        "is_missing_idats": "Missing IDATs",
        "is_missing_gtc": "Missing GTCs",
        "is_cr1_filtered": "Call Rate 1 Filtered",
        "is_cr2_filtered": "Call Rate 2 Filtered",
    }

    if remove_contam:
        exclusion_criteria["is_contaminated"] = "Contamination"

    if remove_rep_discordant:
        exclusion_criteria["is_discordant_replicate"] = "Replicate Discordance"

    sample_qc["analytic_exclusion"] = sample_qc.reindex(exclusion_criteria.keys(), axis=1).any(
        axis=1
    )
    sample_qc["num_analytic_exclusion"] = (
        sample_qc.reindex(exclusion_criteria.keys(), axis=1).sum(axis=1).astype(int)
    )
    sample_qc["analytic_exclusion_reason"] = _get_reason(sample_qc, exclusion_criteria)

    return sample_qc


def _add_subject_representative(sample_qc: pd.DataFrame) -> pd.DataFrame:
    """Flag indicating which sample to use as subject representative.

    First we remove the following samples:

    - Samples to be kept for subject level analysis (``analytic_exclusion``)
    - Internal QC controls (``is_internal_control``)

    For subject IDs with multiple remaining samples, we select the sample
    that has the highest Call Rate 2.
    """
    sub_repr = []
    for _, dd in sample_qc.query("not analytic_exclusion & not is_internal_control").groupby(
        "Group_By_Subject_ID"
    ):
        best_cr = dd.Call_Rate_2.argmax()
        sub_repr.append(dd.iloc[best_cr].name)

    sample_qc["is_subject_representative"] = sample_qc.index.isin(sub_repr)

    return sample_qc


def _add_subject_dropped_from_study(sample_qc: pd.DataFrame) -> pd.Series:
    """Flag indicating which subjects have no representative sample.

    This flag excludes internal controls which by nature are ignored in
    subject level analysis.
    """
    subject_w_no_sample = []  # sourcery skip: list-comprehension
    for sub, dd in sample_qc.query("not is_internal_control").groupby("Group_By_Subject_ID"):
        if dd.is_subject_representative.sum() == 0:
            subject_w_no_sample.append(sub)
    sample_qc["subject_dropped_from_study"] = sample_qc.Group_By_Subject_ID.isin(
        subject_w_no_sample
    )
    return sample_qc


def save(sample_qc: pd.DataFrame, filename: Path) -> None:
    """Save main QC table."""
    sample_qc.reindex(DTYPES, axis=1).to_csv(filename, index=False)


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {"contam": None, "intensity": None}
        defaults.update({k: (Path(v) if v else None) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({k: v for k, v in snakemake.params.items()})  # type: ignore # noqa
        defaults.update({"outfile": Path(snakemake.output[0])})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
