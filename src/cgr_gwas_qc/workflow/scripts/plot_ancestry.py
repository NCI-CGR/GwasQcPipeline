"""
plot_ancestry.py
----------------

This script plots the triangle plot of ancestries from GRAF ancestry
estimates.

Output:
    ``sample_level/ancestry.png``

"""
import os
from pathlib import Path
from typing import Optional

import pandas as pd
import seaborn as sns
import ternary
import typer

from cgr_gwas_qc.reporting import CASE_CONTROL_COLORS
from cgr_gwas_qc.workflow.scripts import sample_qc_table

app = typer.Typer(add_completion=False)


@app.command()
def main(sample_qc: Path, outfile: Path):
    sample = load_sample_data(sample_qc)
    plot(sample, outfile)


def load_sample_data(sample_qc: Path) -> pd.DataFrame:
    return (
        sample_qc_table.read(sample_qc)
        .query("is_subject_representative")
        .dropna(subset=["EUR", "AFR", "ASN"])
    )


def plot(sample: pd.DataFrame, outfile: Optional[os.PathLike] = None):
    sns.set_context("paper")  # use seaborn's context to make sane plot defaults for a paper

    # Create plots
    style_defaults = dict(linewidth=0, alpha=0.8, s=8)
    fig, tax = ternary.figure(scale=1)  # Set scale 0 to 1
    fig.set_size_inches(6, 5)

    # Plot cases, controls, QC, and unknowns separately. Make sure case is last so most visible
    case = sample.query("case_control == 'Case'")
    if case.shape[0] > 0:
        case_color = CASE_CONTROL_COLORS[0]
        tax.scatter(
            case[["EUR", "AFR", "ASN"]].values,
            color=case_color,
            label="Case",
            zorder=4,
            **style_defaults
        )

    control = sample.query("case_control == 'Control'")
    if control.shape[0] > 0:
        control_color = CASE_CONTROL_COLORS[1]  # blue
        tax.scatter(
            control[["EUR", "AFR", "ASN"]].values,
            color=control_color,
            label="Control",
            **style_defaults
        )

    # Issue 237: Add samples if they are neither case or control.
    project_qc = sample.query("case_control == 'QC'")
    if project_qc.shape[0] > 0:
        project_qc_color = CASE_CONTROL_COLORS[2]  # Yellow
        tax.scatter(
            project_qc[["EUR", "AFR", "ASN"]].values,
            color=project_qc_color,
            label="QC",
            **style_defaults
        )

    unknown = sample.query(
        "case_control != 'Control' and case_control != 'Case' and case_control != 'QC'"
    )
    if unknown.shape[0] > 0:
        unknown_color = CASE_CONTROL_COLORS[3]  # Gray
        tax.scatter(
            unknown[["EUR", "AFR", "ASN"]].values,
            color=unknown_color,
            label="Unknown",
            **style_defaults
        )

    # Add plot elements
    multiple = 0.1  # Our scale is 0 to 1 and we want 0.1 increments
    tax.boundary(linewidth=0.5)
    tax.gridlines(multiple=multiple, color="gray")
    tax.ticks(axis="lbr", linewidth=1, multiple=multiple, offset=0.02, tick_formats="%.1f")

    # Set Axis labels
    label_defaults = dict(fontsize=12, offset=0.14)
    tax.left_axis_label("Asian", **label_defaults)
    tax.right_axis_label("African", **label_defaults)
    tax.bottom_axis_label("European", **label_defaults)

    # Add legend
    tax.legend(title="case_control")

    # Clean-up plot
    tax.set_background_color(color="white")
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis("off")  # removes outer square axes

    # Save if given an outfile
    if outfile:
        tax.savefig(outfile)


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({"sample_qc": Path(snakemake.input[0])})  # type: ignore # noqa
        defaults.update({"outfile": Path(snakemake.output[0])})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
