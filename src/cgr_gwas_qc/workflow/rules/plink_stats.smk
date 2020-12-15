"""Plink Stats modules

Calculate a variety of statistics that are supported by plink.

- sexcheck
- allele frequency
- HWE
"""


rule plink_stats_sexcheck:
    input:
        bed="{prefix}/samples.bed",
        bim="{prefix}/samples.bim",
        fam="{prefix}/samples.fam",
    params:
        in_prefix="{prefix}/samples",
        out_prefix="{prefix}/samples",
    output:
        sexcheck="{prefix}/samples.sexcheck",
    group:
        "plink_stats"
    resources:
        mem=10000,
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    shell:
        "plink "
        "--bfile {params.in_prefix} "
        "--check-sex "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"


rule plink_stats_allele_freq:
    input:
        bed="{prefix}/samples.bed",
        bim="{prefix}/samples.bim",
        fam="{prefix}/samples.fam",
    params:
        in_prefix="{prefix}/samples",
        out_prefix="{prefix}/samples",
    output:
        sexcheck="{prefix}/samples.frq",
    group:
        "plink_stats"
    resources:
        mem=10000,
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    shell:
        "plink "
        "--bfile {params.in_prefix} "
        "--freq "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"


rule plink_stats_hardy:
    input:
        bed="{prefix}/samples.bed",
        bim="{prefix}/samples.bim",
        fam="{prefix}/samples.fam",
    params:
        in_prefix="{prefix}/samples",
        out_prefix="{prefix}/samples",
    output:
        sexcheck="{prefix}/samples.hwe",
    group:
        "plink_stats"
    resources:
        mem=10000,
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    shell:
        "plink "
        "--bfile {params.in_prefix} "
        "--hardy "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"
