"""
plot_autosomal_heterozygosity.py
--------------------------------

This script plots the distribution of autosomal heterozygosity for a given population.

Output:
    ``population_level/{population}/autosomal_heterozygosity.png``

"""

from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import typer

from cgr_gwas_qc.parsers.plink import read_het
from cgr_gwas_qc.reporting import CASE_CONTROL_COLORS
from cgr_gwas_qc.workflow.scripts import subject_qc_table

app = typer.Typer(add_completion=False)

COLORS = [
    CASE_CONTROL_COLORS[0],
    CASE_CONTROL_COLORS[1],
    CASE_CONTROL_COLORS[3],
]  # Only want Case/Control/Unknown colors b/c I am dropping QC.


@app.command()
def main(qc_table: Path, het: Path, population: str, threshold: float, outfile: Path):
    df = (
        read_het(het)
        .join(subject_qc_table.read(qc_table).set_index("Group_By_Subject_ID"), how="left")
        .reindex(["case_control", "F"], axis=1)
        .apply(_update_categories)
        .sort_values("F")
        .assign(x_label=lambda x: range(1, x.shape[0] + 1))
    )

    g = plot(df, population, threshold)

    if outfile:
        g.savefig(outfile)


def plot(df: pd.DataFrame, population: str, threshold: float):
    sns.set_context("paper")  # use seaborn's context to make sane plot defaults for a paper

    CASE_CONTROL_LABEL_COLORS = {
        "Case": CASE_CONTROL_COLORS[0],
        "Control": CASE_CONTROL_COLORS[1],
        "QC": CASE_CONTROL_COLORS[2],
        "Unknown": CASE_CONTROL_COLORS[3],
    }

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.scatterplot(
        x="x_label",
        y="F",
        data=df,
        hue="case_control",
        palette=CASE_CONTROL_LABEL_COLORS,
        ax=ax,
        alpha=0.8,
        linewidth=0,
        s=10,
    )

    # Add labels
    ax.set_xlabel("Subjects sorted by F")
    ax.set_ylabel("F")
    ax.set_ylim(_get_ylim(df.F, threshold))
    ax.set_title(f"{population} Heterozygosity F Coefficient")

    # Move legend
    plt.legend(loc="upper left")

    # Add threshold lines
    line_defaults = dict(color="k", ls="--")
    ax.axhline(threshold, **line_defaults)
    ax.axhline(-threshold, **line_defaults)

    # Remove outside edges for a cleaner plot
    sns.despine(ax=ax)

    return fig


def _get_ylim(F: pd.Series, threshold: float) -> Tuple[float, float]:
    default_ylim = threshold + 0.15  # Looks nice when threshold +/- 0.1

    # If there are extreme values outside of default_ylim, then make sure they get plotted
    ylim_min = min(-default_ylim, F.min())
    ylim_max = max(default_ylim, F.max())

    return ylim_min, ylim_max


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
