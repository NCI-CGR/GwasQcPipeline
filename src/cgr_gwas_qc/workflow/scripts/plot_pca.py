"""
plot_pca.py
-----------

This script creates a pair-plot for the first 6 PCs for each population. If
an output directory is given then this script will save a png for each
population.

Output:
    ``population_level/pca_png/{population}.png``

"""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import typer

from cgr_gwas_qc.parsers.eigensoft import Eigenvec
from cgr_gwas_qc.reporting import CASE_CONTROL_COLORS
from cgr_gwas_qc.workflow.scripts import subject_qc_table

app = typer.Typer(add_completion=False)

COLORS = [
    CASE_CONTROL_COLORS[0],
    CASE_CONTROL_COLORS[1],
    CASE_CONTROL_COLORS[3],
]  # Only want Case/Control/Unknown colors b/c I am dropping QC.


@app.command()
def main(qc_table: Path, eigenvec: Path, population: str, outfile: Path):
    pca = (
        Eigenvec(eigenvec)
        .components.join(
            subject_qc_table.read(qc_table).set_index("Group_By_Subject_ID"), how="left"
        )
        .reindex(["case_control", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6"], axis=1)
        .apply(_update_categories)
    )

    g = plot(pca, population)

    if outfile:
        g.savefig(outfile)


def plot(df: pd.DataFrame, population: str) -> sns.PairGrid:
    sns.set_context("paper")  # use seaborn's context to make sane plot defaults for a paper

    CASE_CONTROL_LABEL_COLORS = {'Case': CASE_CONTROL_COLORS[0], 'Control': CASE_CONTROL_COLORS[1], 'QC': CASE_CONTROL_COLORS[2], 'Unknown': CASE_CONTROL_COLORS[3]}

    g = sns.PairGrid(df, hue="case_control", palette=CASE_CONTROL_LABEL_COLORS, corner=True)
    g.map_lower(sns.scatterplot, s=10, alpha=0.8, linewidth=0)
    g.map_diag(sns.kdeplot)
    g.add_legend(
        bbox_to_anchor=(0.2, 0.9), loc=2, prop={"size": "large"}, title_fontsize="xx-large"
    )
    plt.text(0.5, 0.5, population, fontsize=24, rotation=-45, transform=g.fig.transFigure)

    return g


def _update_categories(sr: pd.Series):
    """Update categorical data types for nicer plots"""
    if sr.name == "case_control":
        return sr.cat.remove_categories(["QC"])
    return sr


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update(dict(snakemake.params))  # type: ignore # noqa
        defaults.update({"outfile": Path(snakemake.output[0])})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
