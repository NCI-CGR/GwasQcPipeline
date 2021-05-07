from pathlib import Path
from textwrap import dedent

import pandas as pd


wildcard_constraints:
    name="samples|subjects|controls|snps",
    filters=".*",
    cr="1|2",
    maf="0\.\d+",
    ld="0\.\d+",
    deliver_prefix=".*",
    deliver_suffix=".*",


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
        ped=temp("{prefix}.eigenstrat.ped"),
        map_=temp("{prefix}.eigenstrat.map"),
    log:
        "{prefix}.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2")
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
        "--out {wildcards.prefix}.eigenstrat"


rule concordance_table:
    """Parse IBD file and calculate sample concordance.

    Returns:
        DataFrame with ["ID1", "ID2", "PI_HAT", "concordance"].
          Sample concordance is `IBS2 / (IBS0 + IBS1 + IBS2)`.
    """
    input:
        "{prefix}.genome",
    params:
        concordance_threshold=cfg.config.software_params.dup_concordance_cutoff,
        pi_hat_threshold=cfg.config.software_params.pi_hat_threshold,
    output:
        "{prefix}.concordance.csv",
    script:
        "../scripts/concordance_table.py"


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
        "gen": f"{prefix}.eigenstrat.ped",
        "snp": f"{prefix}.eigenstrat.map",
        "ind": f"{prefix}.eigenstrat.ped",
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
    genotypename: {prefix}.eigenstrat.ped
    snpname: {prefix}.eigenstrat.map
    indivname: {prefix}.eigenstrat.ped
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
        ped="{prefix}.eigenstrat.ped",
        map_="{prefix}.eigenstrat.map",
        par="{prefix}.convert.par",
    output:
        gen=temp("{prefix}.gen"),
        snp=temp("{prefix}.snp"),
        ind=temp("{prefix}.ind"),
    envmodules:
        cfg.envmodules("eigensoft"),
    conda:
        cfg.conda("eigensoft")
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
        cfg.conda("eigensoft")
    shell:
        "smartpca -p {input.par}"
