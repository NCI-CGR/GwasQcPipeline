from pathlib import Path

import pandas as pd
import typer

from cgr_gwas_qc.parsers import sample_sheet
from cgr_gwas_qc.reporting import REPORT_NAME_MAPPER
from cgr_gwas_qc.typing import PathLike
from cgr_gwas_qc.workflow.scripts import population_qc_table, sample_concordance, sample_qc_table

app = typer.Typer(add_completion=False)


@app.command()
def main(
    sample_sheet_csv: Path,
    sample_concordance_csv: Path,
    sample_qc_csv: Path,
    population_qc_csv: Path,
    graf: Path,
    outfile: Path,
):

    with pd.ExcelWriter(outfile) as writer:
        _all_qc(sample_sheet_csv, sample_qc_csv).to_excel(writer, sheet_name="ALL_QC")
        _ancestry(sample_qc_csv, graf).to_excel(writer, sample_sheet="ANCESTRY")
        _concordance(sample_qc_csv, sample_concordance_csv).to_excel(
            writer, sample_sheet="SAMPLE_CONCORDANCE"
        )
        _families(population_qc_csv).to_excel(writer, sample_sheet="FAMILIES")
        _pca(population_qc_csv, writer)
        _het(population_qc_csv, writer)


_ALL_QC_COLUMNS = [
    # Core IDs
    "Sample_ID",
    "Group_By_Subject_ID",
    "Project",
    "Case/Control_Status",
    "Internal Control",
    "Number Samples Per Subject",
    "Replicate IDs",
    "Sample Used in Subject Analysis",
    "Sample Excluded from QC",
    "Subject Removed",
    # QC Results
    "Sample Pass QC",
    "IdatsInProjectDir",
    "Low Call Rate",
    "Call_Rate_Initial",
    "Call_Rate_1_filter",
    "Call_Rate_1",
    "Call_Rate_2_filter",
    "Call_Rate_2",
    "Contaminated",
    "IdatIntensity",
    "Contamination_Rate",
    "Sex Discordant",
    "Expected_Sex",
    "Predicted_Sex",
    "ChrX_Inbreed_estimate",
    "Unexpected Replicate",
    "AFR",
    "EUR",
    "ASN",
    "Ancestry",
    "Count_of_QC_Issue",
]


def _all_qc(sample_sheet_csv: PathLike, sample_qc_csv: PathLike) -> pd.DataFrame:
    ss = sample_sheet.read(sample_sheet_csv).rename(REPORT_NAME_MAPPER, axis=1)
    _additional_columns = [x for x in ss.columns if x not in _ALL_QC_COLUMNS]
    return (
        sample_qc_table.read(sample_qc_csv)
        .rename(REPORT_NAME_MAPPER, axis=1)
        .merge(ss, on="Sample_ID", suffixes=["", "_DROP"])
        .filter(regex="^(?!.*_DROP)")
        .reindex(_ALL_QC_COLUMNS + _additional_columns, axis=1)
    )


_ANCESTRY_COLUMNS = [
    "Sample_ID",
    "Case/Control_Status",
    "#SNPs",
    "GD1",
    "GD2",
    "GD3",
    "GD4",
    "F(%)",
    "E(%)",
    "A(%)",
    "African",
    "European",
    "Asian",
    "Mexican",
    "Indian-Pakistani",
    "Ancestry",
]


def _ancestry(sample_qc_csv: PathLike, graf_txt: PathLike) -> pd.DataFrame:
    sample_qc = sample_qc_table.read(sample_qc_csv).reindex(
        ["Sample_ID", "case_control", "Ancestry"], axis=1
    )

    graf = (
        pd.read_csv(graf_txt, sep="\t")
        .rename({"Sample": "Sample_ID"}, axis=1)
        .drop("DS No.", axis=1)
    )

    return (
        graf.merge(sample_qc, on="Sample_ID")
        .rename(REPORT_NAME_MAPPER, axis=1)
        .reindex(_ANCESTRY_COLUMNS, axis=1)
    )


_CONCORDANCE_COLUMNS = [
    "Sample_ID1",
    "Sample_ID2",
    "Subject_ID1",
    "Subject_ID2",
    "Internal Control1",
    "Internal Control2",
    "Case/Control_Status1",
    "Case/Control_Status2",
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


def _concordance(sample_qc_csv: PathLike, sample_concordance_csv: PathLike) -> pd.DataFrame:
    metadata = (
        sample_qc_table.read(sample_qc_csv)
        .set_index("Sample_ID")
        .rename(REPORT_NAME_MAPPER, axis=1)
        .reindex(["Case/Control_Status", "Ancestry"], axis=1)
    )
    return (
        sample_concordance.read(sample_concordance_csv)
        .merge(metadata.rename_axis("Sample_ID1").add_suffix("1"), on="Sample_ID1", how="left")
        .merge(metadata.rename_axis("Sample_ID2").add_suffix("2"), on="Sample_ID2", how="left")
        .rename(REPORT_NAME_MAPPER, axis=1)
        .rename(
            {
                "is_internal_control1": "Internal Control1",
                "is_internal_control2": "Internal Control2",
            },
            axis=1,
        )
        .reindex(_CONCORDANCE_COLUMNS, axis=1)
    )


_FAMILY_COLUMNS = ["Fam_ID", "Subject_ID"]


def _families(population_qc_csv: PathLike) -> pd.DataFrame:
    return (
        population_qc_table.read(population_qc_csv)
        .dropna(subset=["QC_Family_ID"])
        .rename(REPORT_NAME_MAPPER, axis=1)
        .reindex(_FAMILY_COLUMNS, axis=1)
    )


_PCA_COLUMNS = [
    "Subject_ID",
    "CaCo",
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
    df = population_qc_table.read(population_qc_csv)
    for population, dd in df.groupby("population"):
        dd.reindex(_PCA_COLUMNS, axis=1).to_excel(writer, sheet_name=f"{population}_PCA")


_HET_COLUMNS = ["Subject_ID", "CaCo", "O(HOM)", "E(HOM)", "N(NM)", "F"]


def _het(population_qc_csv: PathLike, writer: pd.ExcelWriter):
    df = population_qc_table.read(population_qc_csv)
    for population, dd in df.groupby("population"):
        dd.reindex(_HET_COLUMNS, axis=1).to_excel(writer, sheet_name=f"{population}_HET")


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {
            **{k: Path(v) for k, v in snakemake.input.items()},  # type: ignore # noqa
            **{"outfile": Path(snakemake.output[0])},  # type: ignore # noqa
        }
        main(**defaults)
    else:
        app()
