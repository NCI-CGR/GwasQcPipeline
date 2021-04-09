from pathlib import Path

import typer

from cgr_gwas_qc.parsers import sample_sheet
from cgr_gwas_qc.reporting import REPORT_NAME_MAPPER
from cgr_gwas_qc.workflow.scripts import sample_qc_table

app = typer.Typer(add_completion=False)

COLUMNS = (
    "SR_Subject_ID",
    "LIMS_Individual_ID",
    "Sample_ID",
    "Project-Sample ID",
    "Call Rate",
    "Low Call Rate",
    "Contaminated",
    "Sex Discordant",
    "Expected Replicate Discordance",
    "Unexpected Replicate",
)


@app.command()
def main(sample_sheet_csv: Path, sample_qc_csv: Path, outfile: Path):
    ss = sample_sheet.read(sample_sheet_csv, all_user_column=True, remove_exclusions=False)
    qc = sample_qc_table.read_sample_qc(sample_qc_csv).rename(REPORT_NAME_MAPPER, axis=1)

    # Merge and drop duplicate rows from ss
    df = qc.merge(ss, on="Sample_ID", how="outer", suffixes=["", "_DROP"]).filter(
        regex="^(?!.*_DROP)", axis=1
    )

    # Adjust names and column order to match legacy
    df.reindex(COLUMNS, axis=1).to_csv(outfile, index=False)


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {
            **{k: Path(v) for k, v in snakemake.input.items()},  # type: ignore # noqa
            **{"outfile": Path(snakemake.output[0])},  # type: ignore # noqa
        }
        main(**defaults)
    else:
        app()
