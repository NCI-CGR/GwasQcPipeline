"""
plot_hwe.py
-----------

This script plots the p-value QQ-plot for the Hardy-Weinberg equilibrium exact-test.

Output:
    ``population_level/hwe_plots/{population}.png``
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import typer

from cgr_gwas_qc.parsers import plink

app = typer.Typer(add_completion=False)


@app.command()
def main(hwe_file: Path, population: str, outfile: Path):
    df = load_p_values(hwe_file)
    g = plot(df, population)

    if outfile:
        g.savefig(outfile)


def load_p_values(filename: Path) -> pd.DataFrame:
    observed_p = plink.read_hwe(filename).P.sort_values()
    expected_p = sorted(
        np.random.uniform(size=observed_p.shape[0])
    )  # p-values should be distributed uniformly
    return pd.DataFrame({"observed_p": observed_p, "expected_p": expected_p})


def plot(df: pd.DataFrame, population: str):
    sns.set_context("paper")  # use seaborn's context to make sane plot defaults for a paper

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(df.expected_p, df.observed_p, color="r")  # add line graph

    # Add labels
    ax.set_xlabel("Expected P value")
    ax.set_ylabel("Observed P value")
    ax.set_title(f"{population} Hardy-Weinberg Proportion")

    # Add diagonal line
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="gray")

    # Add summary counts table at different p-value cutoffs
    ax.text(0.1, 0.1, get_counts(df), transform=ax.transAxes, fontfamily="monospace")

    # Remove outside edges for a cleaner plot
    sns.despine(ax=ax)

    return fig


def get_counts(df: pd.DataFrame) -> str:
    return (
        pd.Series(
            {
                "total": df.shape[0],
                "p ≤ 0.05": (df.observed_p <= 0.05).sum(),
                "p ≤ 0.001": (df.observed_p <= 0.001).sum(),
            }
        )
        .apply("{:,}".format)
        .to_frame()
        .T.to_string(index=False, col_space=15)
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({"hwe_file": Path(snakemake.input[0])})  # type: ignore # noqa
        defaults.update(dict(snakemake.params))  # type: ignore # noqa
        defaults.update({"outfile": Path(snakemake.output[0])})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
