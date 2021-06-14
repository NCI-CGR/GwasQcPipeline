from dataclasses import dataclass
from typing import Dict, Union

import pandas as pd

from cgr_gwas_qc.reporting.templating import number_formater
from cgr_gwas_qc.typing import PathLike
from cgr_gwas_qc.workflow.scripts import population_qc_table, sample_qc_table, subject_qc_table

DataFrameOrFile = Union[pd.DataFrame, PathLike]


@dataclass
class ExclusionTables:
    sample_exclusions: str
    subject_exclusions: str

    @classmethod
    def construct(
        cls,
        sample_qc: DataFrameOrFile,
        subject_qc: DataFrameOrFile,
        population_qc: DataFrameOrFile,
    ) -> "ExclusionTables":

        sample_exclusions = sample_exclusion_counts(sample_qc)
        subject_exclusions = subject_exclusion_counts(subject_qc, population_qc)

        return cls(
            sample_exclusions=to_pretty_markdown(sample_exclusions),
            subject_exclusions=to_pretty_markdown(subject_exclusions),
        )


def sample_exclusion_counts(data: DataFrameOrFile) -> pd.DataFrame:
    if isinstance(data, pd.DataFrame):
        df = data.copy()
    else:
        df = sample_qc_table.read(data)

    return _create_sample_exclusion_counts_table(df)


def _create_sample_exclusion_counts_table(df: pd.DataFrame) -> pd.DataFrame:
    row_names = {
        "is_user_exclusion": "User Excluded",
        "is_array_processing_failure": "Array Processing Failure",  # added below
        "is_cr1_filtered": "Completion Rate 1st Filter",
        "is_cr2_filtered": "Completion Rate 2nd Filter",
        "contam_pass_cr": "Contaminated",
        "internal_control_pass": "Internal QC Samples Removed",
        "samples_remaining": "Samples Remaining for Analysis",  # added below
        "dropped_replicate": "Expected Duplicates Removed",  # added below
        "is_subject_representative": "Subjects Remaining",
    }

    col_names = {
        "Control": "Control",
        "Case": "Case",
        "QC": "QC",
        "Unknown": "Other",
        "total": "All Samples",  # added below
    }

    return (
        df.assign(is_array_processing_failure=lambda x: x.is_missing_idats | x.is_missing_gtc)
        .assign(samples_remaining=lambda x: ~x.analytic_exclusion & ~x.is_internal_control)
        .assign(dropped_replicate=lambda x: x.samples_remaining & ~x.is_subject_representative)
        .assign(contam_pass_cr=lambda x: x.is_contaminated & ~x.is_call_rate_filtered)
        .assign(
            internal_qc_pass=lambda x: x.is_internal_control
            & ~x.is_sample_exclusion
            & ~x.is_contaminated
            & ~x.is_call_rate_filtered
        )
        .reindex(["case_control", *row_names.keys()], axis=1)  # pull out desired columns
        .groupby("case_control")
        .sum()  # Count the number True
        .pipe(_clean_up_counts_table, row_names=row_names, col_names=col_names)
    )


def subject_exclusion_counts(subject: DataFrameOrFile, population: DataFrameOrFile) -> pd.DataFrame:
    # Subject Level (Sex Discordance, Unexpected Replicates)
    subject_df = (
        subject.copy() if isinstance(subject, pd.DataFrame) else subject_qc_table.read(subject)
    )
    subject_exclusions = _create_subject_exclusion_counts_table(subject_df)

    # Population Level (Autosomal Heterozygosity)
    population_df = (
        population.copy()
        if isinstance(population, pd.DataFrame)
        else population_qc_table.read(population)
    )
    population_exclusions = _create_population_exclusion_counts_table(population_df)

    return subject_exclusions.append(population_exclusions).fillna(0)


def _create_subject_exclusion_counts_table(df: pd.DataFrame) -> pd.DataFrame:
    row_names = {
        "is_sex_discordant": "Sex Discordant",
        "is_unexpected_replicate": "Unexpected Replicates",
    }

    col_names = {
        "Control": "Control",
        "Case": "Case",
        "QC": "QC",
        "Unknown": "Other",
        "total": "All Subjects",
    }

    return (
        df.reindex(["case_control", *row_names.keys()], axis=1)
        .groupby("case_control")
        .sum()  # Count the number True
        .pipe(_clean_up_counts_table, row_names=row_names, col_names=col_names)
    )


def _create_population_exclusion_counts_table(df: pd.DataFrame) -> pd.DataFrame:
    row_names = {
        "is_extreme_autosomal_heterozygosity": "Autosomal Het",
    }

    col_names = {
        "Control": "Control",
        "Case": "Case",
        "QC": "QC",
        "Unknown": "Other",
        "total": "All Subjects",
    }

    return (
        df.reindex(["case_control", *row_names.keys()], axis=1)
        .groupby("case_control")
        .sum()
        .pipe(_clean_up_counts_table, row_names=row_names, col_names=col_names)
    )


def _clean_up_counts_table(
    df: pd.DataFrame, row_names: Dict[str, str], col_names: Dict[str, str]
) -> pd.DataFrame:
    return (
        df.reset_index()
        .assign(
            case_control=lambda x: x.case_control.astype("object")
        )  # convert categorical to string: needed for adding the total column
        .set_index("case_control")
        .T.pipe(add_totals, axis=1)  # Add a total column
        .reindex(row_names.keys())  # Reorder and Rename rows and columns
        .rename(row_names)
        .reindex(col_names.keys(), axis=1)
        .rename(col_names, axis=1)
        .rename_axis("Filter Reason/Description")  # re-label the index
    )


def add_totals(df, axis=0, name="total"):
    return df.assign(**{name: lambda x: x.sum(axis).astype(int)})


def to_pretty_markdown(df: pd.DataFrame) -> str:
    return df.applymap(number_formater).to_markdown()
