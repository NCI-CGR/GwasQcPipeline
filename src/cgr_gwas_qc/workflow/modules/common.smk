from pathlib import Path
from textwrap import dedent

import pandas as pd


wildcard_constraints:
    name="samples|subjects|controls|snps",
    filters=".*",
    cr="1|2",
    maf="0\.\d+",
    ld="0\.\d+",


def rule_group(wildcards):
    prefix = wildcards.get("prefix", "").replace("/", "_")
    name = wildcards.get("name", "")
    return f"{prefix}_{name}"


rule plink_bed_to_ped:
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    output:
        ped=temp("{prefix}.ped"),
        map_=temp("{prefix}.map"),
    log:
        "{prefix}.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--recode "
        "--keep-allele-order "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {wildcards.prefix}"


rule concordance_table:
    """Parse IBD file and calculate sample concordance.

    Returns:
        DataFrame with ["IID1", "IID2", "PI_HAT", "concordance"].
          Sample concordance is `IBS2 / (IBS0 + IBS1 + IBS2)`.
    """
    input:
        "{prefix}.genome",
    output:
        temp("{prefix}.concordance.csv"),
    run:
        (
            pd.read_csv(input[0], delim_whitespace=True)
            .assign(concordance=lambda x: x.IBS2 / (x.IBS0 + x.IBS1 + x.IBS2))
            .reindex(["IID1", "IID2", "PI_HAT", "concordance"], axis=1)
            .to_csv(output[0])
        )


def eigensoft_config_inputs(wildcards):
    tool = wildcards.tool
    prefix = wildcards.prefix
    if tool == "pca":
        return {
            "gen": f"{prefix}.gen",
            "snp": f"{prefix}.snp",
            "ind": f"{prefix}.ind",
        }

    # Default is to convert PED/MAP to EIGENSTRAT
    return {
        "gen": f"{prefix}.ped",
        "snp": f"{prefix}.map",
        "ind": f"{prefix}.ped",
    }


def eigensoft_config_params(wildcards):
    tool = wildcards.tool
    prefix = wildcards.prefix

    if tool == "pca":
        return f"""\
        genotypename: {prefix}.gen
        snpname: {prefix}.snp
        indivname: {prefix}.ind
        evecoutname: {prefix}.eigenvec
        fastmode: YES
        """

    # Default is to convert PED/MAP to EIGENSTRAT
    return f"""\
    genotypename: {prefix}.ped
    snpname: {prefix}.map
    indivname: {prefix}.ped
    outputformat: EIGENSTRAT
    genooutfilename: {prefix}.gen
    snpoutfilename: {prefix}.snp
    indoutfilename: {prefix}.ind
    familynames: NO
    """


rule eigensoft_config:
    input:
        unpack(eigensoft_config_inputs),
    params:
        eigensoft_config_params,
    output:
        par=temp("{prefix}.{tool}.par"),
    wildcard_constraints:
        tool="convert|pca",
    run:
        Path(output.par).write_text(dedent(params[0]))


rule eigensoft_convert:
    input:
        ped="{prefix}.ped",
        map_="{prefix}.map",
        par="{prefix}.convert.par",
    output:
        gen=temp("{prefix}.gen"),
        snp=temp("{prefix}.snp"),
        ind=temp("{prefix}.ind"),
    envmodules:
        cfg.envmodules("eigensoft"),
    conda:
        cfg.conda("eigensoft.yml")
    shell:
        "convertf -p {input.par}"


rule eigensoft_smartpca:
    input:
        gen="{prefix}.gen",
        snp="{prefix}.snp",
        ind="{prefix}.ind",
        par="{prefix}.pca.par",
    output:
        gen="{prefix}.eigenvec",
    envmodules:
        cfg.envmodules("eigensoft"),
    conda:
        cfg.conda("eigensoft.yml")
    shell:
        "smartpca -p {input.par}"
