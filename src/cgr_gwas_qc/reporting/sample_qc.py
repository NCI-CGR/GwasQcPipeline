from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from cgr_gwas_qc.reporting.templating import number_formater


@dataclass
class SampleQC:
    array_processing: "ArrayProcessing"
    completion_rate: "CompletionRate"
    contamination: "Contamination"
    internal_controls: "InternalControls"
    expected_replicates: "ExpectedReplicates"
    summary: "SampleSummary"
    ancestry: "Ancestry"

    @classmethod
    def construct(
        cls,
        sample_sheet: pd.DataFrame,
        snp_qc: pd.DataFrame,
        sample_qc: pd.DataFrame,
        control_replicates: pd.DataFrame,
        study_replicates: pd.DataFrame,
        call_rate_png: Path,
        ancestry_png: Path,
    ) -> "SampleQC":
        return cls(
            ArrayProcessing.construct(sample_sheet),
            CompletionRate.construct(snp_qc, sample_qc, call_rate_png),
            Contamination.construct(sample_qc),
            InternalControls.construct(sample_qc, control_replicates),
            ExpectedReplicates.construct(sample_qc, study_replicates),
            SampleSummary.construct(sample_qc),
            Ancestry.construct(sample_qc, ancestry_png),
        )


@dataclass
class ArrayProcessing:
    num_samples_excluded: int
    num_missing_idats: int
    num_missing_gtc: int
    num_array_processing_failure: int
    num_user_exclusions: int
    num_samples_qc_processed: int

    @classmethod
    def construct(cls, sample_sheet: pd.DataFrame) -> "ArrayProcessing":
        return cls(
            num_samples_excluded=sample_sheet.is_sample_exclusion.sum(),
            num_missing_idats=sample_sheet.is_missing_idats.sum(),
            num_missing_gtc=sample_sheet.is_missing_gtc.sum(),
            num_array_processing_failure=(
                sample_sheet.is_missing_idats | sample_sheet.is_missing_gtc
            ).sum(),
            num_user_exclusions=sample_sheet.is_user_exclusion.sum(),
            num_samples_qc_processed=(~sample_sheet.is_sample_exclusion).sum(),
        )


@dataclass
class CompletionRate:
    mean_initial_call_rate: float

    num_snp_cr1_filtered: int
    num_snp_pass_cr1: int

    num_sample_cr1_filtered: int
    num_sample_pass_cr1: int

    num_snp_cr2_filtered: int
    num_snp_pass_cr2: int

    num_sample_cr2_filtered: int
    num_sample_pass_cr2: int

    num_sample_pass_call_rate: int
    num_snp_pass_call_rate: int

    png: str

    @classmethod
    def construct(
        cls, snp_qc: pd.DataFrame, sample_qc: pd.DataFrame, png: Path
    ) -> "CompletionRate":
        return cls(
            mean_initial_call_rate=sample_qc.Call_Rate_Initial.mean(),
            num_snp_cr1_filtered=snp_qc.is_cr1_filtered.sum(),
            num_snp_pass_cr1=(~snp_qc.is_cr1_filtered).sum(),
            num_sample_cr1_filtered=sample_qc.is_cr1_filtered.sum(),
            num_sample_pass_cr1=(~sample_qc.is_cr1_filtered).sum(),
            num_snp_cr2_filtered=snp_qc.is_cr2_filtered.sum(),
            num_snp_pass_cr2=(~snp_qc.is_cr2_filtered).sum(),
            num_sample_cr2_filtered=sample_qc.is_cr2_filtered.sum(),
            num_sample_pass_cr2=(~sample_qc.is_cr2_filtered).sum(),
            num_sample_pass_call_rate=(~sample_qc.is_call_rate_filtered).sum(),
            num_snp_pass_call_rate=(~snp_qc.is_call_rate_filtered).sum(),
            png=png.resolve().as_posix(),
        )


@dataclass
class Contamination:
    num_pass_cr_and_contaminated: int
    num_remaining: int

    @classmethod
    def construct(cls, sample_qc: pd.DataFrame) -> "Contamination":
        return cls(
            sample_qc.query("not is_call_rate_filtered & is_contaminated").shape[0],
            sample_qc.query("not is_call_rate_filtered & not is_contaminated").shape[0],
        )


@dataclass
class InternalControls:
    num_internal_controls: int
    num_internal_controls_pass: int
    num_remaining: int
    min_concordance: float
    mean_concordance: float

    @classmethod
    def construct(cls, sample_qc: pd.DataFrame, replicates: pd.DataFrame):
        # List of internal controls that passed QC steps.
        ic_pass = sample_qc.query(  # type: ignore # noqa
            "is_internal_control"
            " & not is_sample_exclusion"
            " & not is_call_rate_filtered"
            " & not is_contaminated"
        ).Sample_ID.tolist()

        return cls(
            sample_qc.is_internal_control.sum(),
            sample_qc.query("Sample_ID == @ic_pass").shape[0],
            sample_qc.query(
                "not is_sample_exclusion"
                " & not is_call_rate_filtered"
                " & not is_contaminated"
                " & not is_internal_control"
            ).shape[0],
            replicates.query(
                "Sample_ID1 == @ic_pass & Sample_ID2 == @ic_pass"
            ).PLINK_concordance.min(),
            replicates.query(
                "Sample_ID1 == @ic_pass & Sample_ID2 == @ic_pass"
            ).PLINK_concordance.mean(),
        )


@dataclass
class ExpectedReplicates:
    num_low_concordance: int
    num_not_subject_representative: int
    num_remaining: int
    min_concordance: float
    mean_concordance: float

    @classmethod
    def construct(cls, sample_qc: pd.DataFrame, replicates: pd.DataFrame) -> "ExpectedReplicates":
        return cls(
            sample_qc["Expected Replicate Discordance"].sum(),
            sample_qc.query(
                "not is_internal_control & not analytic_exclusion & not is_subject_representative"
            ).shape[0],
            sample_qc.query("is_subject_representative").shape[0],
            replicates.PLINK_concordance.min(),
            replicates.PLINK_concordance.mean(),
        )


@dataclass
class SampleSummary:
    case_control: str
    male_female: str
    num_starting_subjects: int
    num_samples_excluded: int
    num_subjects_remaining: int

    @classmethod
    def construct(cls, sample_qc: pd.DataFrame) -> "SampleSummary":
        return cls(
            case_control="{}:{}".format(
                *sample_qc.case_control.value_counts()[["Case", "Control"]]
            ),
            male_female="{}:{}".format(*sample_qc.expected_sex.value_counts()[["M", "F"]]),
            num_starting_subjects=cls._count_subjects(sample_qc),
            num_samples_excluded=sample_qc.query("not is_subject_representative").shape[0],
            num_subjects_remaining=sample_qc.query(
                "not is_internal_control & is_subject_representative"
            ).shape[0],
        )

    @staticmethod
    def _count_subjects(sample_qc: pd.DataFrame) -> int:
        Subject_IDs = sample_qc.Group_By_Subject_ID
        Subject_IDs[sample_qc.is_internal_control] = (
            "QC"  # Treat internal controls as a single subject
        )
        return Subject_IDs.unique().shape[0]


@dataclass
class Ancestry:
    table: str
    png: str

    @classmethod
    def construct(cls, sample_qc: pd.DataFrame, png: Path) -> "Ancestry":
        return cls(cls._build_table(sample_qc), png.resolve().as_posix())

    @staticmethod
    def _build_table(sample_qc: pd.DataFrame) -> str:
        # Issue 281: Ensure 'Other' ancestry output from graf-pop captured in table 3

        # Add 'Other' to the categories and full blanks with "Other"
        sample_qc["Ancestry"] = sample_qc["Ancestry"].fillna("Other")

        return (
            sample_qc.query("is_subject_representative")
            .groupby("Ancestry")
            .size()
            .rename("Count")
            .rename_axis("Ancestry Group")
            .map(number_formater)
            .sort_index()
            .to_markdown()
        )
