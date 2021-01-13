"""Plink Stats modules

Calculate a variety of statistics that are supported by plink.

- sexcheck
- allele frequency
- HWE
"""


rule plink_stats_call_rate:
    """Runs ``plink`` missingness statistics.

    .. rubric:: Statistics
        :imiss`: sample-based missing data report
        :lmiss`: variant-based missing data report
    """
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    output:
        imiss="{prefix}.imiss",
        lmiss="{prefix}.lmiss",
    group:
        "plink_stats"
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 2
    resources:
        mem=10000,
    shell:
        "plink "
        "--bfile {wildcards.prefix} "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--missing "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {wildcards.prefix}"


rule plink_stats_sexcheck:
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    output:
        sexcheck="{prefix}.sexcheck",
    group:
        "plink_stats"
    threads: 2
    resources:
        mem=10000,
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--check-sex "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {wildcards.prefix}"


rule plink_stats_allele_freq:
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    output:
        sexcheck="{prefix}.frq",
    group:
        "plink_stats"
    threads: 2
    resources:
        mem=10000,
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--freq "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {wildcards.prefix}"


rule plink_stats_hardy:
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    output:
        sexcheck="{prefix}.hwe",
    group:
        "plink_stats"
    threads: 2
    resources:
        mem=10000,
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--hardy "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {wildcards.prefix}"
