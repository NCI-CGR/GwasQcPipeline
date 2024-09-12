#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""Extracts Illumina scores from a multisample VCF for target sample and writes into a single sample adpc.bin format."""
from pathlib import Path
from typing import Generator

import numpy as np
import pandas as pd
import typer
from pysam import VariantFile

# from more_itertools import unzip
from cgr_gwas_qc.parsers.illumina import AdpcRecord, AdpcWriter

# In[ ]:


app = typer.Typer(add_completion=False)


@app.command()
def main(
    bcf_file: Path = typer.Argument(
        ..., help="Path to a multisample VCF file with needed scores.", exists=True, readable=True
    ),
    target_sample: str = typer.Argument(..., help="Sample name."),
    outfile: Path = typer.Argument(
        ..., help="Path to output the adpc.bin file.", file_okay=True, writable=True
    ),
) -> None:
    """Extracts Illumina scores from a multisample VCF for target sample and writes into a single sample adpc.bin format.

    The adpc.bin format includes::

        x_raw: Raw intensities for allele A
        y_raw: Raw intensities for allele B
        x_norm: Normalized intensities for allele A
        y_norm: Normalized intensities for allele B
        genotype_score: The genotype clustering confidence score
        genotype: The called genotype (0: AA, 1: AB, 2: BB, 3: unknown or missing)
    """

    vcf = VariantFile(bcf_file)
    sample_fields_to_fetch = ["X", "Y", "NORMX", "NORMY", "GT"]
    info_fieilds_to_fetch = ["ALLELE_A", "ALLELE_B", "GC_SCORE"]

    renaming_scheme = {
        "X": "x_raw",
        "Y": "y_raw",
        "NORMX": "x_norm",
        "NORMY": "y_norm",
        "GC_SCORE": "genotype_score",
    }

    vcf_info = []
    for rec in vcf:
        vcf_info.append(
            [rec.info.get(key) for key in info_fieilds_to_fetch]
            + [rec.samples[target_sample].get(key) for key in sample_fields_to_fetch]
        )
    vcf_info = pd.DataFrame(vcf_info, columns=info_fieilds_to_fetch + sample_fields_to_fetch)

    # convert GT to illumina genotype:
    vcf_info = pd.concat(
        [
            vcf_info.reset_index(drop=True),
            pd.DataFrame(vcf_info["GT"].tolist(), columns=["GT1", "GT2"]),
        ],
        axis=1,
    )
    vcf_info["genotype"] = np.nan  # setting missing by default.
    vcf_info.loc[
        (vcf_info["GT1"] == vcf_info["ALLELE_A"]) & (vcf_info["GT2"] == vcf_info["ALLELE_A"]),
        "genotype",
    ] = "AA"
    vcf_info.loc[
        (vcf_info["GT1"] == vcf_info["ALLELE_B"]) & (vcf_info["GT2"] == vcf_info["ALLELE_B"]),
        "genotype",
    ] = "BB"
    vcf_info.loc[(vcf_info["GT1"] != vcf_info["GT2"]), "genotype"] = "AB"
    vcf_info.loc[(vcf_info["GT1"].isna()) | (vcf_info["GT2"].isna()), "genotype"] = np.nan

    # renaming for consisting with gtc2adpc
    vcf_info.rename(columns=renaming_scheme, inplace=True)

    # encode illumina genotypes per adpc.bin encoding
    vcf_info = vcf_info.replace(
        {"genotype": {"AA": 0, "AB": 1, "BA": 1, "BB": 2, "A": 3, "B": 3}}
    ).fillna({"genotype": 3})
    vcf_info["genotype_score"] = vcf_info["genotype_score"].where(vcf_info["x_raw"] != 0, 0)

    # writing adpc.bin file
    with AdpcWriter(outfile) as fh:
        for record in vcf_to_adpc_record(vcf_info):
            fh.write(record)


# In[2]:


## function to adpc records
def vcf_to_adpc_record(vcf_info) -> Generator[AdpcRecord, None, None]:
    for x_raw, y_raw, x_norm, y_norm, genotype_score, genotype in zip(
        vcf_info["x_raw"],
        vcf_info["y_raw"],
        vcf_info["x_norm"],
        vcf_info["y_norm"],
        vcf_info["genotype_score"],
        vcf_info["genotype"],
    ):
        yield AdpcRecord(x_raw, y_raw, x_norm, y_norm, genotype_score, genotype)


# In[ ]:


# def adpc_as_df(adpc):
#     i=0
#     adpc_records = []
#     with AdpcReader(adpc) as fh:
#                    for record in fh:
#                       # if i <10000:
#                            adpc_records.append(record)
#                            i+=1
#     return(pd.DataFrame(adpc_records,columns=['x_raw', 'y_raw', 'x_norm', 'y_norm', 'genotype_score', 'genotype']))


# In[ ]:


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({k: v for k, v in snakemake.params.items()})  # type: ignore # noqa
        defaults.update({"outfile": Path(snakemake.output[0])})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
