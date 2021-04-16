from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from cgr_gwas_qc.models.config import Config
from cgr_gwas_qc.reporting.templating import number_formater


@dataclass
class SampleQC:
    array_processing: "ArrayProcessing"
    completion_rate: "CompletionRate"
    contamination: "Contamination"
    internal_controls: "InternalControls"
    expected_replicates: "ExpectedReplicates"
    unexpected_replicates: "UnExpectedReplicates"
    summary: "SampleSummary"
    sex_verification: "SexVerification"
    ancestry: "Ancestry"

    @classmethod
    def construct(
        cls,
        config: Config,
        sample_sheet: pd.DataFrame,
        snp_qc: pd.DataFrame,
        sample_qc: pd.DataFrame,
        control_replicates: pd.DataFrame,
        study_replicates: pd.DataFrame,
        unexpected_replicates: pd.DataFrame,
        call_rate_png: Path,
        chrx_inbreeding_png: Path,
        ancestry_png: Path,
    ) -> "SampleQC":

        return cls(
            ArrayProcessing.construct(config, sample_sheet),
            CompletionRate.construct(snp_qc, sample_qc, call_rate_png),
            Contamination.construct(sample_qc),
            InternalControls.construct(sample_qc, control_replicates),
            ExpectedReplicates.construct(sample_qc, study_replicates),
            UnExpectedReplicates.construct(sample_qc, unexpected_replicates),
            SampleSummary.construct(sample_qc),
            SexVerification.construct(sample_qc, chrx_inbreeding_png),
            Ancestry.construct(sample_qc, ancestry_png),
        )


@dataclass
class ArrayProcessing:
    num_failed_array_processing: int
    num_samples_qc_processed: int

    @classmethod
    def construct(cls, config: Config, sample_sheet: pd.DataFrame) -> "ArrayProcessing":
        if config.Sample_IDs_to_remove:
            num_failed_array_processing = len(config.Sample_IDs_to_remove)
        else:
            num_failed_array_processing = 0

        return cls(
            num_failed_array_processing=num_failed_array_processing,
            num_samples_qc_processed=sample_sheet.shape[0],
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
    num_remaining: int
    min_concordance: float
    mean_concordance: float

    @classmethod
    def construct(cls, sample_qc: pd.DataFrame, replicates: pd.DataFrame):
        return cls(
            sample_qc.is_internal_control.sum(),
            sample_qc.query(
                "not is_call_rate_filtered " "& not is_contaminated " "& not is_internal_control"
            ).shape[0],
            replicates.PLINK_concordance.min(),
            replicates.PLINK_concordance.mean(),
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
            sample_qc.is_replicate_discordant.sum(),
            sample_qc.query("is_pass_sample_qc & not is_subject_representative").shape[0],
            sample_qc.query("is_subject_representative").shape[0],
            replicates.PLINK_concordance.min(),
            replicates.PLINK_concordance.mean(),
        )


@dataclass
class UnExpectedReplicates:
    num_unexpected_replicates: int
    min_concordance: float
    mean_concordance: float

    @classmethod
    def construct(cls, sample_qc: pd.DataFrame, replicates: pd.DataFrame) -> "UnExpectedReplicates":
        return cls(
            sample_qc.is_unexpected_replicate.sum(),
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
        Subject_IDs[
            sample_qc.is_internal_control
        ] = "QC"  # Treat internal controls as a single subject
        return Subject_IDs.unique().shape[0]


@dataclass
class SexVerification:
    num_sex_discordant: int
    num_remaining: int
    png: str

    @classmethod
    def construct(cls, sample_qc: pd.DataFrame, png: Path) -> "SexVerification":
        return cls(
            sample_qc.query("is_subject_representative & is_sex_discordant").shape[0],
            sample_qc.query("is_subject_representative & not is_sex_discordant").shape[0],
            png.resolve().as_posix(),
        )


@dataclass
class Ancestry:

    table: str
    png: str

    @classmethod
    def construct(cls, sample_qc: pd.DataFrame, png: Path) -> "Ancestry":
        return cls(cls._build_table(sample_qc), png.resolve().as_posix())

    @staticmethod
    def _build_table(sample_qc: pd.DataFrame) -> str:
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
