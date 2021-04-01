#!/usr/bin/env python
import os
from pathlib import Path

import numpy as np
import pandas as pd
import typer

from cgr_gwas_qc import load_config
from cgr_gwas_qc.config import ConfigMgr
from cgr_gwas_qc.paths import make_path_list
from cgr_gwas_qc.reporting import env
from cgr_gwas_qc.validators import check_file

app = typer.Typer(add_completion=False)

template = env.get_template("summary_stats.txt.j2")


@app.command()
def main(file_name: Path, file_out: Path):
    """Generate summary statistics from the sample qc table.

    Args:
        file_name (Path): Path to the `all_sample_qc.csv` file.
        file_out (Path): Path to the output file.
    """
    cfg = load_config()
    df = pd.read_csv(file_name)

    payload = {
        "working_dir": os.getcwd(),
        "n_samples": df.shape[0],
        "n_no_gtc": _number_missing_gtc_files(cfg, df),
        "idats": _value_counts_w_na(df.idats_exist),
        "call_rate_initial": _table_summary_str(df.Call_Rate_Initial),
        "contamination_rate": _table_summary_w_na(df.Contamination_Rate),
        "intensity": _table_summary_str(df.IdatIntensity),
        "contaminated": _value_counts_w_na(df.is_contaminated),
        "call_rate_filtered": _value_counts_w_na(df["is_call_rate_filtered"]),
        "sex_discord": _value_counts_w_na(df["is_sex_discordant"]),
        "rep_discord": _value_counts_w_na(df["is_replicate_discordant"]),
        "unexpected_rep": _value_counts_w_na(df["is_unexpected_replicate"]),
        "identifiler": _value_counts_w_na(df["identifiler_needed"]),
        "issues": _value_counts_w_na(df["Count_of_QC_Issue"]),
    }

    file_out.write_text(template.render(**payload))


def _number_missing_gtc_files(cfg: ConfigMgr, df: pd.DataFrame) -> int:
    if cfg.config.user_files.gtc_pattern is None:
        return df.shape[0]

    cnt = 0
    for file_name in make_path_list(cfg.expand(cfg.config.user_files.gtc_pattern)):
        try:
            check_file(file_name)
        except (FileNotFoundError, PermissionError):
            cnt += 1

    return cnt


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
