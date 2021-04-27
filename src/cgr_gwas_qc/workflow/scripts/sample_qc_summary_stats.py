#!/usr/bin/env python
import os
from pathlib import Path

import numpy as np
import pandas as pd
import typer

from cgr_gwas_qc.reporting import env
from cgr_gwas_qc.workflow.scripts import sample_qc_table

app = typer.Typer(add_completion=False)

template = env.get_template("summary_stats.txt.j2")


@app.command()
def main(filename: Path, outfile: Path):
    """Generate summary statistics from the sample qc table.
    """
    df = sample_qc_table.read(filename)

    payload = {
        "working_dir": os.getcwd(),
        "n_samples": df.shape[0],
        "n_no_gtc": df.is_missing_gtc.sum(),
        "idats": _value_counts_w_na(df.is_missing_idats),
        "call_rate_initial": _table_summary_str(df.Call_Rate_Initial),
        "contamination_rate": _table_summary_w_na(df.Contamination_Rate),
        "intensity": _table_summary_str(df.IdatIntensity),
        "contaminated": _value_counts_w_na(df.is_contaminated),
        "call_rate_filtered": _value_counts_w_na(df["is_call_rate_filtered"]),
        "rep_discord": _value_counts_w_na(df["is_discordant_replicate"]),
        "identifiler": _value_counts_w_na(df["identifiler_needed"]),
        "issues": _value_counts_w_na(df["num_analytic_exclusion"]),
    }

    outfile.write_text(template.render(**payload))


def _table_summary(sr: pd.Series) -> pd.DataFrame:
    return (
        sr.describe()
        .map(lambda x: np.round(x, 4))
        .to_frame()
        .T[["min", "25%", "50%", "mean", "75%", "max"]]
    )


def _table_summary_str(sr: pd.Series) -> str:
    return _table_summary(sr).to_string(index=False)


def _table_summary_w_na(sr: pd.Series) -> str:
    return _table_summary(sr).assign(NaN=sr.isna().sum()).to_string(index=False)


def _value_counts_w_na(sr: pd.Series) -> str:
    return (
        sr.value_counts(dropna=False).sort_index().to_frame().T.to_string(na_rep="NA", index=False)
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        main(Path(snakemake.input[0]), Path(snakemake.output[0]))  # type: ignore # noqa
    else:
        app()
