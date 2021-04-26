from dataclasses import dataclass

import pandas as pd

from cgr_gwas_qc.models.config import Config
from cgr_gwas_qc.reporting.templating import number_formater


@dataclass
class ExclusionTables:
    sample_exclusions: str
    subject_exclusions: str

    @classmethod
    def construct(
        cls,
        config: Config,
        sample_sheet: pd.DataFrame,
        sample_qc: pd.DataFrame,
        population_qc: pd.DataFrame,
    ) -> "ExclusionTables":

        pre_qc = cls._pre_qc(config, sample_sheet)
        pop_qc = cls._pop_qc(population_qc)
        df = (
            sample_qc.merge(pre_qc, how="outer")
            .merge(pop_qc, how="left")
            .merge(sample_qc, how="left")
        )

        return cls(
            sample_exclusions=cls._sample_exclusion_table(df),
            subject_exclusions=cls._subject_exclusion_table(df),
        )

    @staticmethod
    def _pre_qc(config: Config, sample_sheet: pd.DataFrame) -> pd.DataFrame:
        df = sample_sheet.copy()
        df["is_array_processing_failure"] = False

        if config.Sample_IDs_to_remove:
            mask = df.Sample_ID.isin(config.Sample_IDs_to_remove)
            df.loc[mask, "is_array_processing_failure"] = True

        return df.reindex(["Sample_ID", "is_array_processing_failure"], axis=1)

    @staticmethod
    def _pop_qc(population_qc: pd.DataFrame) -> pd.DataFrame:
        return population_qc.reindex(
            ["Subject_ID", "is_extreme_autosomal_heterozygosity"], axis=1
        ).rename({"Subject_ID": "Group_By_Subject_ID"}, axis=1)

    @staticmethod
    def _sample_exclusion_table(df: pd.DataFrame) -> str:
        row_names = {
            "is_array_processing_failure": "Array Processing Failure",
            "is_cr1_filtered": "Completion Rate 1st Filter",
            "is_cr2_filtered": "Completion Rate 2nd Filter",
            "is_contaminated": "Contaminated",
            "is_internal_control": "Internal QC Samples Removed",
            "samples_remaining": "Samples Remaining for Analysis",
            "dropped_replicate": "Expected Duplicates Removed",  # added below
            "is_subject_representative": "Subjects Remaining",
        }

        col_names = {
            "Control": "Controls",
            "Case": "Cases",
            "QC": "QC",
            "Unknown": "Other",
            "total": "All Samples",
        }

        agg_df = (
            df.assign(samples_remaining=lambda x: ~x.analytic_exclusion & ~x.is_internal_control)
            .assign(dropped_replicate=lambda x: x.samples_remaining & ~x.is_subject_representative)
            .reindex(["case_control", *row_names.keys()], axis=1)
            .groupby("case_control")
            .sum()
        )
        agg_df.index = agg_df.index.astype("object")  # remove categorical

        return (
            agg_df.T.assign(total=lambda x: x.sum(axis=1).astype("Int32"))
            .reindex(row_names.keys())
            .rename(row_names)
            .reindex(col_names.keys(), axis=1)
            .rename(col_names, axis=1)
            .rename_axis("Filter Reason/Description")
            .applymap(number_formater)
            .to_markdown()
        )

    @staticmethod
    def _subject_exclusion_table(df: pd.DataFrame) -> str:
        row_names = {
            "is_sex_discordant": "Sex Discordant",
            "is_unexpected_replicate": "Unexpected Replicates",
            "is_extreme_autosomal_heterozygosity": "Autosomal Het",
        }

        col_names = {
            "Control": "Controls",
            "Case": "Cases",
            "QC": "QC",
            "Unknown": "Other",
            "total": "All Subjects",
        }

        agg_df = (
            df.query("is_subject_representative")
            .reindex(["case_control", *row_names.keys()], axis=1)
            .groupby("case_control")
            .sum()
        )
        agg_df.index = agg_df.index.astype("object")  # remove categorical

        return (
            agg_df.T.assign(total=lambda x: x.sum(axis=1).astype("Int32"))
            .reindex(row_names.keys())
            .rename(row_names)
            .reindex(col_names.keys(), axis=1)
            .rename(col_names, axis=1)
            .rename_axis("Filter Reason")
            .applymap(number_formater)
            .to_markdown()
        )
