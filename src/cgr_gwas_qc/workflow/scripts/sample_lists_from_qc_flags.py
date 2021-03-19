#!/usr/bin/env python
from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False)


def main(
    all_samples: Path = typer.Argument(..., help="Sample level QC report."),
    cr: Path = typer.Argument(..., help="Path to output file low_call_rate.txt"),
    contam: Path = typer.Argument(..., help="Path to output file contaminated.txt"),
    sex: Path = typer.Argument(..., help="Path to output file sex_discordant.txt"),
    rep: Path = typer.Argument(..., help="Path to output file replicate_discordant.txt"),
    ctrl: Path = typer.Argument(..., help="Path to output file internal_controls.txt"),
):

    df = pd.read_csv(all_samples)
    _save_sample_flag_as_file(df, "call_rate_filtered", cr)
    _save_sample_flag_as_file(df, "Contaminated", contam)
    _save_sample_flag_as_file(df, "sex_discordant", sex)
    _save_sample_flag_as_file(df, "Expected Replicate Discordance", rep)
    _save_sample_flag_as_file(df, "internal_control", ctrl)


def _save_sample_flag_as_file(df, col, file_name):
    (
        df.fillna(False)  # If missing assume False
        .assign(Sample_ID2=lambda x: x.Sample_ID)
        .query(f"`{col}`")
        .reindex(["Sample_ID", "Sample_ID2"], axis=1)
        .to_csv(file_name, index=False, header=False, sep=" ")
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({k: Path(v) for k, v in snakemake.output.items()})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
