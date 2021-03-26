"""
plot_pca.py
-----------

This script creates a pair-plot for the first 6 PCs for each population. If
an output directory is given then this script will save a png for each
population.

Output:
    ``population_level/pca_png/{population}.png``

"""
import os
from pathlib import Path
from typing import Optional

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
def main(population_qc: Path, outdir: Path):
    outdir.mkdir(exist_ok=True, parents=True)

    df = load_population_data(population_qc)
    for population, dd in df.groupby("population"):
        outfile = outdir / f"{population}.png"
        plot(dd, population, outfile)


def load_population_data(filename: Path) -> pd.DataFrame:
    return (
        read_agg_population_qc_tables(filename)
        .query("case_control != 'QC'")
        .transform(_update_categories)
        .reindex(["population", "case_control", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6"], axis=1)
    )


def plot(df: pd.DataFrame, population: str, outfile: Optional[os.PathLike] = None):
    sns.set_context("paper")  # use seaborn's context to make sane plot defaults for a paper

    g = sns.PairGrid(df, hue="case_control", palette=COLORS, corner=True)
    g.map_lower(sns.scatterplot, s=10, alpha=0.8, linewidth=0)
    g.map_diag(sns.kdeplot)
    g.add_legend(
        bbox_to_anchor=(0.2, 0.9), loc=2, prop={"size": "large"}, title_fontsize="xx-large"
    )
    plt.text(0.5, 0.5, population, fontsize=24, rotation=-45, transform=g.fig.transFigure)

    if outfile:
        g.savefig(outfile)


def _update_categories(sr: pd.Series):
    """Update categorical data types for nicer plots"""
    if sr.name == "case_control":
        return sr.cat.remove_categories(["QC"])
    return sr


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({"population_qc": Path(snakemake.input[0])})  # type: ignore # noqa
        defaults.update({"outdir": Path(snakemake.output[0])})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
