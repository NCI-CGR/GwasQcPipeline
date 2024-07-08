from pathlib import Path

import pandas as pd
import typer

from cgr_gwas_qc.parsers import sample_sheet
from cgr_gwas_qc.reporting import REPORT_NAME_MAPPER
from cgr_gwas_qc.typing import PathLike
from cgr_gwas_qc.workflow.scripts import (
    agg_population_concordance,
    population_qc_table,
    sample_concordance,
    sample_qc_table,
    subject_qc_table,
)

app = typer.Typer(add_completion=False)


@app.command()
def main(
    sample_sheet_csv: Path,
    sample_concordance_csv: Path,
    sample_qc_csv: Path,
    subject_qc_csv: Path,
    population_concordance_csv: Path,
    population_qc_csv: Path,
    graf: Path,
    outfile: Path,
):

    with pd.ExcelWriter(outfile) as writer:
        _sample_qc(sample_sheet_csv, sample_qc_csv).to_excel(
            writer, sheet_name="SAMPLE_QC", index=False
        )
        _sample_concordance(sample_qc_csv, sample_concordance_csv).to_excel(
            writer, sheet_name="SAMPLE_CONCORDANCE", index=False
        )
        _subject_qc(sample_sheet_csv, subject_qc_csv).to_excel(
            writer, sheet_name="SUBJECT_QC", index=False
        )
        _ancestry(sample_qc_csv, graf).to_excel(writer, sheet_name="ANCESTRY", index=False)
        _families(population_qc_csv, writer)
        _population_concordance(population_concordance_csv, writer)
        _pca(population_qc_csv, writer)
        _het(population_qc_csv, writer)


_SAMPLE_QC_COLUMNS = [
    # Core IDs
    "Sample_ID",
    "Group_By_Subject_ID",
    "Project",
    "Case/Control_Status",
    "Internal Control",
    "Number Samples Per Subject",
    "Replicate IDs",
    "Sample Excluded from QC",
    "IdatsInProjectDir",
    "Sample Used in Subject Analysis",
    "Subject Removed",
    # QC Results
    "Sample Pass QC",
    "Low Call Rate",
    "Call_Rate_Initial",
    "Call_Rate_1_filter",
    "Call_Rate_1",
    "Call_Rate_2_filter",
    "Call_Rate_2",
    "Contaminated",
    "Contamination_Rate",
    "IdatIntensity",
    "Expected Replicate",
    "Expected Replicate Discordance",
    "Unexpected Replicate",
    "Sex Discordant",
    "Expected_Sex",
    "Predicted_Sex",
    "ChrX_Inbreed_estimate",
    "AFR",
    "EUR",
    "ASN",
    "Ancestry",
    "Count_of_QC_Issue",
]


def _sample_qc(sample_sheet_csv: PathLike, sample_qc_csv: PathLike) -> pd.DataFrame:
    ss = sample_sheet.read(sample_sheet_csv, remove_exclusions=False).rename(
        REPORT_NAME_MAPPER, axis=1
    )
    _additional_columns = [x for x in ss.columns if x not in _SAMPLE_QC_COLUMNS]
    return (
        sample_qc_table.read(sample_qc_csv)
        .rename(REPORT_NAME_MAPPER, axis=1)
        .merge(ss, on="Sample_ID", suffixes=["", "_DROP"], how="outer")
        .filter(regex="^(?!.*_DROP)")
        .reindex(_SAMPLE_QC_COLUMNS + _additional_columns, axis=1)
    )


_SAMPLE_CONCORDANCE_COLUMNS = [
    "Sample_ID1",
    "Sample_ID2",
    "Subject_ID1",
    "Subject_ID2",
    "Case/Control_Status1",
    "Case/Control_Status2",
    "Sample Used in Subject Analysis1",
    "Sample Used in Subject Analysis2",
    "Ancestry1",
    "Ancestry2",
    "Expected Replicate",
    "Expected Replicate Discordance",
    "Unexpected Replicate",
    "PLINK_PI_HAT",
    "PLINK_concordance",
    "PLINK_is_ge_pi_hat",
    "PLINK_is_ge_concordance",
    "GRAF_HGMR",
    "GRAF_AGMR",
    "GRAF_relationship",
    "KING_Kinship",
    "KING_relationship",
]


def _excel_limit_filter(df: pd.DataFrame) -> pd.DataFrame:
    """Remove extra rows if larger than the excel limit."""
    max_excel_rows = 1_000_000  # the real max is 1_048_576
    if df.shape[0] < max_excel_rows:
        return df

    return df.sort_values("PLINK_PI_HAT", ascending=False).head(max_excel_rows)


def _sample_concordance(sample_qc_csv: PathLike, sample_concordance_csv: PathLike) -> pd.DataFrame:
    ancestry = sample_qc_table.read(sample_qc_csv).set_index("Sample_ID").Ancestry
    representative = (
        sample_qc_table.read(sample_qc_csv).set_index("Sample_ID").is_subject_representative
    )

    return (
        sample_concordance.read(sample_concordance_csv)
        .pipe(
            lambda x: x[
                x.PLINK_PI_HAT.notna() | (x.GRAF_relationship.notna() & x.KING_relationship.notna())
            ]
        )  # Keep if we have results from plink or both GRAF and KING. This will limit the number of rows.
        .merge(ancestry.rename_axis("Sample_ID1").rename("Ancestry1"), on="Sample_ID1", how="left")
        .merge(ancestry.rename_axis("Sample_ID2").rename("Ancestry2"), on="Sample_ID2", how="left")
        .merge(
            representative.rename_axis("Sample_ID1").rename("Sample Used in Subject Analysis1"),
            on="Sample_ID1",
            how="left",
        )
        .merge(
            representative.rename_axis("Sample_ID2").rename("Sample Used in Subject Analysis2"),
            on="Sample_ID2",
            how="left",
        )
        .rename(REPORT_NAME_MAPPER, axis=1)
        .rename(
            {
                "case_control1": REPORT_NAME_MAPPER["case_control"] + "1",
                "case_control2": REPORT_NAME_MAPPER["case_control"] + "2",
            },
            axis=1,
        )
        .reindex(_SAMPLE_CONCORDANCE_COLUMNS, axis=1)
        .pipe(_excel_limit_filter)
        .sort_values(["Sample_ID1", "Sample_ID2"])
    )


_SUBJECT_QC_COLUMNS = [
    "Group_By_Subject_ID",
    "Sample_ID",
    "Case/Control_Status",
    "Unexpected Replicate",
    "unexpected_replicate_ids",
    "Expected_Sex",
    "Predicted_Sex",
    "ChrX_Inbreed_estimate",
    "Sex Discordant",
    "AFR",
    "EUR",
    "ASN",
    "Ancestry",
]


def _subject_qc(sample_sheet_csv: PathLike, subject_qc_csv: PathLike) -> pd.DataFrame:
    ss = sample_sheet.read(sample_sheet_csv).rename(REPORT_NAME_MAPPER, axis=1)
    _additional_columns = [x for x in ss.columns if x not in _SUBJECT_QC_COLUMNS]

    return (
        subject_qc_table.read(subject_qc_csv)
        .rename(REPORT_NAME_MAPPER, axis=1)
        .merge(ss, on="Sample_ID", suffixes=["", "_DROP"])
        .filter(regex="^(?!.*_DROP)")
        .reindex(_SUBJECT_QC_COLUMNS + _additional_columns, axis=1)
    )


_ANCESTRY_COLUMNS = [
    "Sample_ID",
    "Group_By_Subject_ID",
    "Sample Used in Subject Analysis",
    "Case/Control_Status",
    "#SNPs",
    "GD1 (x)",
    "GD2 (y)",
    "GD3 (z)",
    "GD4",
    "F(%)",
    "E(%)",
    "A(%)",
]


def _ancestry(sample_qc_csv: PathLike, graf_txt: PathLike) -> pd.DataFrame:
    sample_qc = sample_qc_table.read(sample_qc_csv).reindex(
        [
            "Sample_ID",
            "Group_By_Subject_ID",
            "case_control",
            "is_subject_representative",
            "Ancestry",
        ],
        axis=1,
    )

    graf = pd.read_csv(graf_txt, sep="\t", skiprows=7).rename({"Sample": "Sample_ID"}, axis=1)

    return (
        graf.merge(sample_qc, on="Sample_ID")
        .rename(REPORT_NAME_MAPPER, axis=1)
        .reindex(_ANCESTRY_COLUMNS, axis=1)
    )


_FAMILY_COLUMNS = ["Fam_ID", "Subject_ID"]


def _families(population_qc_csv: PathLike, writer):
    df = (
        population_qc_table.read(population_qc_csv)
        .dropna(subset=["QC_Family_ID"])
        .rename(REPORT_NAME_MAPPER, axis=1)
    )

    if df.shape[0] > 0:
        for population, dd in df.groupby("population"):
            dd = dd.reindex(_FAMILY_COLUMNS, axis=1)
            dd = dd.sort_values(by=["Fam_ID"])
            dd.to_excel(writer, sheet_name=f"{population[:24]}_FAMILY", index=False)


_POPULATION_CONCORDANCE_COLUMNS = [
    "Subject_ID1",
    "Subject_ID2",
    "Case/Control_Status1",
    "Case/Control_Status2",
    "PLINK_PI_HAT",
    "PLINK_is_ge_pi_hat",
]


def _population_concordance(population_concordance_csv: PathLike, writer: pd.ExcelWriter):
    df = agg_population_concordance.read(population_concordance_csv).rename(
        {
            "case_control1": REPORT_NAME_MAPPER["case_control"] + "1",
            "case_control2": REPORT_NAME_MAPPER["case_control"] + "2",
        },
        axis=1,
    )

    for population, dd in df.groupby("population"):
        dd = (
            dd.reindex(_POPULATION_CONCORDANCE_COLUMNS, axis=1)
            .pipe(_excel_limit_filter)
            .sort_values(["Subject_ID1", "Subject_ID2"])
        )
        dd = dd.sort_values(by=["PLINK_PI_HAT"], ascending=False)
        dd.to_excel(writer, sheet_name=f"{population[:24]}_IBD", index=False)


_PCA_COLUMNS = [
    "Subject_ID",
    "Case/Control_Status",
    "PC1",
    "PC2",
    "PC3",
    "PC4",
    "PC5",
    "PC6",
    "PC7",
    "PC8",
    "PC9",
    "PC10",
]


def _pca(population_qc_csv: PathLike, writer: pd.ExcelWriter):
    df = population_qc_table.read(population_qc_csv).rename(REPORT_NAME_MAPPER, axis=1)
    for population, dd in df.groupby("population"):
        dd.reindex(_PCA_COLUMNS, axis=1).to_excel(
            writer, sheet_name=f"{population[:24]}_PCA", index=False
        )


_HET_COLUMNS = ["Subject_ID", "Case/Control_Status", "O_HOM", "E_HOM", "N_NM", "F"]


def _het(population_qc_csv: PathLike, writer: pd.ExcelWriter):
    df = population_qc_table.read(population_qc_csv).rename(REPORT_NAME_MAPPER, axis=1)
    for population, dd in df.groupby("population"):
        dd.reindex(_HET_COLUMNS, axis=1).to_excel(
            writer, sheet_name=f"{population[:24]}_HET", index=False
        )


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {
            **{k: Path(v) for k, v in snakemake.input.items()},  # type: ignore # noqa
            **{"outfile": Path(snakemake.output[0])},  # type: ignore # noqa
        }
        main(**defaults)
    else:
        app()
