"""outfile: Optional[os.PathLike]
plot_chrx_inbreeding.py
-----------------------

This script plots the distribution of X chromosome inbreeding coefficients by
sex.

Output:
    ``sample_level/chrx_inbreeding.png``

"""
import os
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import typer

# import snakemake

from cgr_gwas_qc.reporting import CASE_CONTROL_COLORS
from cgr_gwas_qc.workflow.scripts import sample_qc_table

app = typer.Typer(add_completion=False)


@app.command()
def main(sample_qc: Path, outfile: Path, xchr: str):
    sample = load_sample_data(sample_qc)
    xchr = str(snakemake.params)  # type: ignore # noqa
    plot(sample, outfile, xchr)


"""
def load_sample_data(sample_qc: Path) -> pd.DataFrame:
    return (
        sample_qc_table.read(sample_qc)
        .query("expected_sex != 'U' and is_subject_representative")  # Don't plot unknown sex
        .transform(_update_categories)
    )
"""


def load_sample_data(sample_qc: Path) -> pd.DataFrame:
    return (
        sample_qc_table.read(sample_qc)
        .query("is_subject_representative")
        .transform(_update_categories)
    )


def _update_categories(sr: pd.DataFrame):
    """Update categorical data types for nicer plots"""
    if sr.name == "case_control":
        print("sr.name == case_control")
        return sr.cat.remove_unused_categories()

    if sr.name == "expected_sex":
        # Drop the 'U' category and re-order to put females first.
        print("sr.name == expected_sex")
        temp = sr.cat.remove_unused_categories()
        print(temp.unique())
        return sr.cat.remove_unused_categories()

    return sr


def plot(sample: pd.DataFrame, outfile: Optional[os.PathLike] = None, xchr: bool = True):
    sns.set_context("paper")  # use seaborn's context to make sane plot defaults for a paper

    # Create plots
    style_defaults = dict(linewidth=0, alpha=0.8, s=2)
    defaults = dict(x="expected_sex", y="X_inbreeding_coefficient", data=sample)
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.boxplot(ax=ax, showfliers=False, **defaults)
    sns.stripplot(
        ax=ax, hue="case_control", palette=CASE_CONTROL_COLORS, **defaults, **style_defaults
    )

    # Make boxplot black and white
    plt.setp(ax.artists, edgecolor="k", facecolor="w")
    plt.setp(ax.lines, color="k")

    # Rename Axes
    # ax.set_xlabel("Reported Sex")
    ax.set_ylabel("ChrX Inbreeding Coeff")

    xchr = xchr.strip().lower() == "true"
    print(type(xchr), " ", xchr)
    if xchr:
        print("sex chr included", xchr)
        ax.set_xlabel("Reported Sex")
    else:
        print("No sex chromosome ", xchr)
        ax.set_xlabel("No sex chromosome \nSkipping sex condordace")

    # Add line at 0.5
    line_defaults = dict(color="k", ls="--")
    ax.axhline(0.5, **line_defaults)

    # Remove lines around inner plots
    sns.despine(ax=ax)

    # Save if given an outfile
    if outfile:
        fig.savefig(outfile)


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({"sample_qc": Path(snakemake.input[0])})  # type: ignore # noqa
        defaults.update({"outfile": Path(snakemake.output[0])})  # type: ignore # noqa
        defaults.update({"xchr": snakemake.params})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
