"""
plot_autosomal_heterozygosity.py
--------------------------------

This script plots the distribution of autosomal heterozygosity for a given population.

Output:
    ``population_level/{population}/autosomal_heterozygosity.png``

"""
import os
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import typer

from cgr_gwas_qc.reporting import CASE_CONTROL_COLORS
from cgr_gwas_qc.workflow.scripts.agg_population_qc_tables import read_agg_population_qc_tables

app = typer.Typer(add_completion=False)

COLORS = [
    CASE_CONTROL_COLORS[0],
    CASE_CONTROL_COLORS[1],
    CASE_CONTROL_COLORS[3],
]  # Only want Case/Control/Unknown colors b/c I am dropping QC.


@app.command()
def main(population_qc: Path, threshold: float, outdir: Path):
    outdir.mkdir(exist_ok=True, parents=True)

    df = load_population_data(population_qc)

    ylim = df.F.agg(["min", "max"]).tolist()
    for population, dd in df.groupby("population"):
        outfile = outdir / f"{population}.png"
        dd = add_x_label(dd)
        plot(dd, population, threshold, ylim, outfile)


def load_population_data(filename: Path) -> pd.DataFrame:
    return (
        read_agg_population_qc_tables(filename)
        .query("case_control != 'QC'")
        .transform(_update_categories)
        .reindex(["population", "case_control", "F"], axis=1)
    )


def add_x_label(df: pd.DataFrame) -> pd.DataFrame:
    """Add X value.

    I want to plot subjects ordered by ``F``. Here I sort by ``F`` and label
    the Subjects from 1 to N.
    """
    return df.sort_values("F").assign(x_label=lambda x: range(1, x.shape[0] + 1))


def plot(
    df: pd.DataFrame,
    population: str,
    threshold: float,
    ylim: Sequence[float],
    outfile: Optional[os.PathLike] = None,
):
    sns.set_context("paper")  # use seaborn's context to make sane plot defaults for a paper

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.scatterplot(
        x="x_label",
        y="F",
        data=df,
        hue="case_control",
        palette=COLORS,
        ax=ax,
        alpha=0.8,
        linewidth=0,
        s=10,
    )

    # Add labels
    ax.set_xlabel("Subjects sorted by F")
    ax.set_ylabel("F")
    ax.set_ylim(ylim)
    ax.set_title(f"{population} Homozygosity F Coefficient")

    # Move legend
    plt.legend(loc="upper left")

    # Add threshold lines
    line_defaults = dict(color="k", ls="--")
    ax.axhline(threshold, **line_defaults)
    ax.axhline(-threshold, **line_defaults)

    # Remove outside edges for a cleaner plot
    sns.despine(ax=ax)

    if outfile:
        fig.savefig(outfile)


def _update_categories(sr: pd.Series):
    """Update categorical data types for nicer plots"""
    if sr.name == "case_control":
        return sr.cat.remove_categories(["QC"])
    return sr


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({"population_qc": Path(snakemake.input[0])})  # type: ignore # noqa
        defaults.update({k: v for k, v in snakemake.params.items()})  # type: ignore # noqa
        defaults.update({"outdir": Path(snakemake.output[0])})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
