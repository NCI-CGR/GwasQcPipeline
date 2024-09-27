"""
Internal SNP QC Report
----------------------

**Script**: ``workflow/scripts/snp_qc_table.py``

This script takes various outputs from the GwasQcPipeline and
aggregates/summarizes results into an internal SNP level QC report.

**Output**: ``sample_level/snp_qc.csv``

.. csv-table::
    :header: name, dtype, description

    array_snp_id, string, The SNP ID as found in the BPM file.
    thousand_genomes_snp_id, string, The rsID found in the 1000 Genomes VCF.
    chromosome, string, The chromosome that the SNP is on.
    snp_cr1, float, The SNP call rate results after the first filter.
    snp_cr1_removed, boolean, True if the SNP was removed by the first call rate filter.
    snp_cr2, float, The SNP call rate results after the second filter.
    snp_cr2_removed, boolean, True if the SNP was removed by the second call rate filter.

"""

from pathlib import Path

import pandas as pd
import typer

from cgr_gwas_qc.parsers import plink

app = typer.Typer(add_completion=False)

DTYPES = {
    "array_snp_id": "string",
    "thousand_genomes_snp_id": "string",
    "chromosome": "string",
    "Call_Rate_Initial": "float",
    "is_cr1_filtered": "boolean",
    "Call_Rate_1": "float",
    "is_cr2_filtered": "boolean",
    "Call_Rate_2": "float",
    "is_call_rate_filtered": "boolean",
}


@app.command()
def main(
    initial: Path = typer.Option(..., help="Initial lmiss file."),
    cr1: Path = typer.Option(..., help="lmiss file after CR1 filters."),
    cr2: Path = typer.Option(..., help="lmiss file after CR2 filters."),
    thousand_genomes: Path = typer.Option(..., help="Mapping of Array IDs to thousand genome IDs."),
    outfile: Path = typer.Option(..., help="Output CSV."),
):

    snp_qc = aggregate_snp_data(initial, cr1, cr2, thousand_genomes)
    add_call_rate_flags(snp_qc)
    snp_qc.reindex(DTYPES.keys(), axis=1).to_csv(outfile, index=False)


def aggregate_snp_data(initial, cr1, cr2, thousand_genomes) -> pd.DataFrame:
    return (
        _read_lmiss(initial, "Call_Rate_Initial")
        .merge(_read_lmiss(cr1, "Call_Rate_1"), how="outer")
        .merge(_read_lmiss(cr2, "Call_Rate_2"), how="outer")
        .merge(_read_thousand_genomes(thousand_genomes), how="outer")
    )


def add_call_rate_flags(df: pd.DataFrame):
    """Add summary flags for call rate.

    Each call rate step will drop snps/samples, so missing snps/samples were
    filtered in the current or previous step(s). For example, a missing
    snp `Call_Rate_1` indicates the snp did not pass Call Rate 1 filter
    or was missing from before. If the sample was missing from before than I
    am setting it as missing (pd.NA). This will help with data provenance and
    hopefully make it clearer why a sample was removed.
    """
    cri = df.Call_Rate_Initial.isna()
    cr1 = df.Call_Rate_1.isna()
    cr2 = df.Call_Rate_2.isna()

    df["is_cr1_filtered"] = cr1
    df.loc[cri, "is_cr1_filtered"] = pd.NA

    df["is_cr2_filtered"] = cr2
    df.loc[cri | cr1, "is_cr2_filtered"] = pd.NA

    # Add Call Rate Summary Flag
    # `True` if filtered in cr1 or cr2
    # `False` if not filtered in cr1 or cr2
    # `pd.NA` if missing in the initial data set (cri)
    df["is_call_rate_filtered"] = cr2
    df.loc[cri, "is_call_rate_filtered"] = pd.NA


def read(filename: Path) -> pd.DataFrame:
    return pd.read_csv(filename, dtype=DTYPES)


def _read_lmiss(filename: Path, name: str) -> pd.DataFrame:
    return (
        plink.read_lmiss(filename)
        .assign(**{name: lambda x: 1 - x.F_MISS})
        .rename({"SNP": "array_snp_id", "CHR": "chromosome"}, axis=1)
        .reindex(["array_snp_id", "chromosome", name], axis=1)
    )


def _read_thousand_genomes(filename: Path):
    return pd.read_csv(filename).rename(
        {"array_id": "array_snp_id", "thousand_genomes_id": "thousand_genomes_snp_id"}, axis=1
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type:ignore # noqa
        defaults.update({"outfile": Path(snakemake.output[0])})  # type:ignore # noqa
        main(**defaults)
    else:
        app()
