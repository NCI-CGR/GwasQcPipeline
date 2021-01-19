from pathlib import Path
from textwrap import dedent

from cgr_gwas_qc.parsers.illumina import BeadPoolManifest


numSNPs = BeadPoolManifest(cfg.config.reference_files.illumina_manifest_file).num_loci


rule plink_bed_to_ped:
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    output:
        ped="{prefix}.ped",
        map_="{prefix}.map",
    log:
        "{prefix}.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 20
    resources:
        mem=10000,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--recode "
        "--keep-allele-order "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {wildcards.prefix}"


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
        par="{prefix}.{tool}.par",
    wildcard_constraints:
        tool="convert|pca",
    group:
        "{tool}"
    run:
        Path(output.par).write_text(dedent(params[0]))


rule eigensoft_convert:
    input:
        ped="{prefix}.ped",
        map_="{prefix}.map",
        par="{prefix}.convert.par",
    output:
        gen="{prefix}.gen",
        snp="{prefix}.snp",
        ind="{prefix}.ind",
    group:
        "convert"
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
    group:
        "convert"
    envmodules:
        cfg.envmodules("eigensoft"),
    conda:
        cfg.conda("eigensoft.yml")
    shell:
        "smartpca -p {input.par}"
