#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import typer
from pysam import VariantFile

app = typer.Typer(add_completion=False)


# In[2]:


def is_biallelic_snp(rec):
    return len(rec.ref) == 1 and len(rec.alts[0]) == 1


# In[27]:


@app.command()
def main(
    vcf_file: Path = typer.Argument(
        ..., help="Path to the sample VCF.", exists=True, readable=True
    ),
    kgvcf_file: Path = typer.Argument(
        ..., help="Path to the 1000 genomes VCF file to get ABF.", exists=True, readable=True
    ),
    population: str = typer.Argument(
        ...,
        help="Which population in VCF file to use for ABF [Not implemented,global/AF used currently].",
    ),
    abf_file: Optional[Path] = typer.Argument(
        None, help="Output file name. If none given it will use the `VCF_FILE + .abf.txt`."
    ),
):
    """Pull out the B allele frequency for a given subpopulation.

    Takes variants from the provided `vcf_file` and tries to find their B
    allele frequency in the provided `kgvcf_file`. The B allele frequencies are
    based on the allele frequency for the provided subpopulation
    `population`. If the provided subpopulation doesn't exist, it will output
    a list of available populations. If the variant does not exists in the
    `vcf_file` then the ABF will be set to "NA".
    """

    if abf_file is None:
        abf_file = vcf_file.with_suffix(".abf.txt")

    # vcf_file='/scratch/rajwanir2/data_catalog/QC-test/results/GSAMD-24v1-0-contamseries/vcf/GSAMD-24v1-0-contamseries.bcf'
    # abf_file='/scratch/rajwanir2/data_catalog/QC-test/GSAMD-24v1-0.new.abf'
    # kg_vcf='/DCEG/CGF/Bioinformatics/Production/data/thousG/hg38/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz'

    # retrieving needed info from the input vcf
    vcf = VariantFile(vcf_file, drop_samples=True)
    vcf_info = []
    for rec in vcf:
        vcf_info.append(
            [
                rec.id,
                rec.chrom,
                rec.pos,
                rec.ref[0],
                rec.alts[0],
                is_biallelic_snp(rec),
                rec.info["GC_SCORE"],
                rec.info["ALLELE_A"],
            ]
        )

    vcf_info = pd.DataFrame(
        vcf_info,
        columns=["id", "chrom", "pos", "ref", "alt", "is_biallelic_snp", "gc_score", "allele_a"],
    )
    vcf_info["is_chrompos_uniq"] = (
        vcf_info.duplicated(subset=["chrom", "pos"], keep=False) == False
    )  # if the input vairants are split, this will catch the multi-allelics

    # fetching B-allele frequencies from the kgvcf
    kg_vcf = VariantFile(kgvcf_file, drop_samples=True)
    kg_info = []
    for i in range(vcf_info.shape[0]):
        if kg_vcf.is_valid_reference_name(vcf_info.chrom[i]):
            for rec in kg_vcf.fetch(
                contig=vcf_info.chrom[i], start=vcf_info.pos[i] - 1, stop=vcf_info.pos[i]
            ):
                if (
                    is_biallelic_snp(rec)
                    & (vcf_info.is_biallelic_snp[i])
                    & (vcf_info.ref[i] == rec.ref[0])
                    & (vcf_info.alt[i] == rec.alts[0])
                    & (vcf_info.is_chrompos_uniq[i])
                ):
                    kg_info.append(
                        [
                            rec.chrom,
                            rec.pos,
                            rec.ref[0],
                            rec.alts[0],
                            rec.info["AC"][0],
                            rec.info["AN"],
                        ]
                    )
        #         else:
        #             kg_info.append([np.nan,np.nan])
        # else:
        #     kg_info.append([np.nan,np.nan])
    kg_info = pd.DataFrame(kg_info, columns=["chrom", "pos", "ref", "alt", "ac", "an"])

    # left-joining the kgvcf B-allele freq to input vcf metrics.
    vcf_info = vcf_info.merge(
        right=kg_info, how="left", sort=False, on=["chrom", "pos", "ref", "alt"]
    )
    # A pseudo way to switch from population AF to Illumina A/B frequency
    vcf_info["af"] = abs(vcf_info["allele_a"] - (vcf_info["ac"] / vcf_info["an"]))
    # Setting AF to na if gentrain score is not >0, if this is egt file was an older one
    vcf_info["af"].where(vcf_info.gc_score > 0, np.nan, inplace=True)
    # writing abf file
    vcf_info[["id", "af"]].to_csv(
        abf_file, sep="\t", index=False, header=["SNP_ID", "ABF"], na_rep="NA"
    )


# In[ ]:


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({k: Path(v) for k, v in snakemake.output.items()})  # type: ignore # noqa
        defaults.update({k: v for k, v in snakemake.params.items()})  # type: ignore # noqa
        main(**defaults)  # type: ignore
    else:
        app()
