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
    "Call_Rate_Initial",
    "Low Call Rate",
    "Contaminated",
    "Sex Discordant",
    "Expected Replicate Discordance",
    "Unexpected Replicate",
)


@app.command()
def main(sample_sheet_csv: Path, sample_qc_csv: Path, outfile: Path):
    ss = sample_sheet.read(sample_sheet_csv, all_user_column=True, remove_exclusions=False)
    qc = sample_qc_table.read(sample_qc_csv).rename(REPORT_NAME_MAPPER, axis=1)

    # Merge and drop duplicate rows from ss
    df = qc.merge(ss, on="Sample_ID", how="outer", suffixes=["", "_DROP"]).filter(
        regex="^(?!.*_DROP)", axis=1
    )
    df = df.reindex(COLUMNS, axis=1)
    # rename columns to work with LIMS system
    df = df.rename(columns = {'Sample_ID':'Sample ID','Call_Rate_Initial':'Call Rate'})
    # For the samples missing GTC files, their status in the low call rate column should be TRUE instead of blank
    df.loc[df['Call Rate'] =="",'Low Call Rate'] = True
    # Fill in blanks for replicate and discordant columns
    df['Unexpected Replicate'] = df['Unexpected Replicate'].replace("",False).fillna(False) 
    df['Expected Replicate Discordance'] = df['Expected Replicate Discordance'].replace("",False).fillna(False)
    df['Sex Discordant'] = df['Sex Discordant'].replace("",False).fillna(False)
    # Adjust names and column order to match legacy
    df.to_csv(outfile, index=False)


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {
            **{k: Path(v) for k, v in snakemake.input.items()},  # type: ignore # noqa
            **{"outfile": Path(snakemake.output[0])},  # type: ignore # noqa
        }
        main(**defaults)
    else:
        app()
