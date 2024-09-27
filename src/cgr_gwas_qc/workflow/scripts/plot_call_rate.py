"""
plot_call_rate.py
-----------------

This script summarizes call rate by plotting SNP and Sample call rate
distributions.

Output:
    ``sample_level/call_rate.png``

"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import typer

from cgr_gwas_qc.reporting import CASE_CONTROL_COLORS
from cgr_gwas_qc.workflow.scripts import sample_qc_table, snp_qc_table

app = typer.Typer(add_completion=False)


@app.command()
def main(
    sample_qc: Path,
    snp_qc: Path,
    sample_cr1: float,
    snp_cr1: float,
    sample_cr2: float,
    snp_cr2: float,
    outfile: Path,
):
    sample = load_sample_data(sample_qc)
    snp = load_snp_data(snp_qc)
    fig = plot_panel(sample, snp, sample_cr1, snp_cr1, sample_cr2, snp_cr2)
    fig.savefig(outfile)


def load_sample_data(sample_qc: Path) -> pd.DataFrame:
    return (
        sample_qc_table.read(sample_qc)
        .reindex(["Call_Rate_Initial", "Call_Rate_2", "case_control"], axis=1)
        .dropna(subset=["Call_Rate_Initial"])
        .transform(_scale_call_rate)
        .sort_values("Call_Rate_Initial")
        .reset_index(drop=True)
        .assign(x_initial=lambda x: x.index / x.shape[0] * 100)
        .sort_values("Call_Rate_2")
        .reset_index(drop=True)
        .assign(x_cr2=lambda x: x.index / x.shape[0] * 100)
    )


def load_snp_data(snp_qc: Path) -> pd.DataFrame:
    return (
        snp_qc_table.read(snp_qc)
        .reindex(["Call_Rate_Initial", "Call_Rate_2"], axis=1)
        .dropna(subset=["Call_Rate_Initial"])
        .transform(_scale_call_rate)
        .sort_values("Call_Rate_Initial")
        .reset_index(drop=True)
        .assign(x_initial=lambda x: x.index / x.shape[0] * 100)
        .sort_values("Call_Rate_2")
        .reset_index(drop=True)
        .assign(x_cr2=lambda x: x.index / x.shape[0] * 100)
    )


def plot_panel(
    sample: pd.DataFrame,
    snp: pd.DataFrame,
    sample_cr1: float,
    snp_cr1: float,
    sample_cr2: float,
    snp_cr2: float,
):
    sns.set_context("paper")  # use seaborn's context to make sane plot defaults for a paper

    fig, axes = plt.subplots(
        2,  # 2 rows of panels
        2,  # 2 columns of panels
        sharex="col",  # share x-axis by columns
        sharey="row",  # share y-axis by rows
        figsize=(6, 8),  # set the size of the figure
        gridspec_kw=dict(hspace=0.1, wspace=0.1),  # set the spacing between figures
    )

    # Set basic defaults so I don't have to repeat myself
    CASE_CONTROL_LABEL_COLORS = {
        "Case": CASE_CONTROL_COLORS[0],
        "Control": CASE_CONTROL_COLORS[1],
        "QC": CASE_CONTROL_COLORS[2],
        "Unknown": CASE_CONTROL_COLORS[3],
    }

    style_defaults = dict(linewidth=0, alpha=0.8, s=5)
    sample_defaults = {
        **dict(hue="case_control", palette=CASE_CONTROL_LABEL_COLORS, data=sample),
        **style_defaults,
    }
    snp_defaults = {**dict(data=snp, palette="gray"), **style_defaults}

    # Add figures Initial Call Rates
    sns.scatterplot(x="x_initial", y="Call_Rate_Initial", ax=axes[0][0], **sample_defaults)
    sns.scatterplot(x="x_initial", y="Call_Rate_Initial", ax=axes[0][1], **snp_defaults)

    # Add figures Call Rates after filtering
    sns.scatterplot(x="x_cr2", y="Call_Rate_2", ax=axes[1][0], **sample_defaults)
    sns.scatterplot(x="x_cr2", y="Call_Rate_2", ax=axes[1][1], **snp_defaults)

    # Add call rate threshold locations
    line_defaults = dict(color="k", ls="--")
    axes[0][0].axhline(sample_cr1 * 100, **line_defaults)
    axes[0][0].axhline(sample_cr2 * 100, **line_defaults)

    axes[0][1].axhline(snp_cr1 * 100, **line_defaults)
    axes[0][1].axhline(snp_cr2 * 100, **line_defaults)

    # Rename X-Axis
    axes[1][0].set_xlabel("Samples (%)")
    axes[1][1].set_xlabel("Loci (%)")

    # Rename Y-Axis
    axes[0][0].set_ylabel("Completion Rate (%)")
    axes[1][0].set_ylabel("Completion Rate (%)\n≥Call Rate 2 Threshold")

    # Set row Y-Axis limits
    axes[0][0].set_ylim(0, None)
    axes[1][0].set_ylim(sample_cr2 * 100, None)

    # Add Title to top row of plots
    axes[0][0].set_title("Sample Completion", fontweight="bold")
    axes[0][1].set_title("Loci Completion", fontweight="bold")

    # Remove lines around inner plots
    sns.despine(ax=axes[0][0])
    sns.despine(ax=axes[0][1], left=True)
    sns.despine(ax=axes[1][0])
    sns.despine(ax=axes[1][1], left=True)

    return fig


def _scale_call_rate(sr: pd.Series) -> pd.Series:
    """Scales call rates from a proportion (0-1) to percentage (0-100)"""
    if sr.name.startswith("Call_Rate"):
        return sr * 100
    return sr


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({"outfile": Path(snakemake.output[0])})  # type: ignore # noqa
        defaults.update({k: float(v) for k, v in snakemake.params.items()})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
