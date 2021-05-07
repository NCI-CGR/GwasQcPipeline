"""Plink Stats modules

Calculate a variety of statistics that are supported by plink.

- sexcheck
- allele frequency
- HWE
"""


include: cfg.modules("common")


rule plink_stats_call_rate:
    """Runs ``plink`` missingness statistics.

    .. rubric:: Statistics
        :imiss`: sample-based missing data report
        :lmiss`: variant-based missing data report
    """
    input:
        bed="{prefix}/{name}{filters}.bed",
        bim="{prefix}/{name}{filters}.bim",
        fam="{prefix}/{name}{filters}.fam",
    params:
        out_prefix="{prefix}/{name}{filters}",
    output:
        imiss="{prefix}/{name}{filters}.imiss",
        lmiss="{prefix}/{name}{filters}.lmiss",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    shell:
        "plink "
        "--bfile {wildcards.prefix} "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--missing "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule plink_stats_sexcheck:
    input:
        bed="{prefix}/{name}{filters}.bed",
        bim="{prefix}/{name}{filters}.bim",
        fam="{prefix}/{name}{filters}.fam",
    params:
        out_prefix="{prefix}/{name}{filters}",
    output:
        sexcheck="{prefix}/{name}{filters}.sexcheck",
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2")
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--check-sex "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule plink_stats_allele_freq:
    input:
        bed="{prefix}/{name}{filters}.bed",
        bim="{prefix}/{name}{filters}.bim",
        fam="{prefix}/{name}{filters}.fam",
    params:
        out_prefix="{prefix}/{name}{filters}",
    output:
        sexcheck="{prefix}/{name}{filters}.frq",
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2")
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--freq "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule plink_stats_hardy:
    input:
        bed="{prefix}/{name}{filters}.bed",
        bim="{prefix}/{name}{filters}.bim",
        fam="{prefix}/{name}{filters}.fam",
    params:
        out_prefix="{prefix}/{name}{filters}",
    output:
        sexcheck="{prefix}/{name}{filters}.hwe",
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2")
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--hardy "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule plink_stats_ibd:
    """Calculates Identity-by-descent.

    This step is trying to find relationships among the samples by
    calculating IBS/IBD. This calculation is no LD-aware so we are doing the
    LD pruning before. This calculation excludes non-autosomes.
    """
    input:
        bed="{prefix}/{name}{filters}.bed",
        bim="{prefix}/{name}{filters}.bim",
        fam="{prefix}/{name}{filters}.fam",
    params:
        ibd_min=cfg.config.software_params.ibd_pi_hat_min,
        ibd_max=cfg.config.software_params.ibd_pi_hat_max,
        out_prefix="{prefix}/{name}{filters}",
    output:
        genome="{prefix}/{name}{filters}.genome",
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2")
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--genome full "
        "--min {params.ibd_min} "
        "--max {params.ibd_max} "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule plink_stats_het:
    """Calculates autosomal heterozygosity."""
    input:
        bed="{prefix}/{name}{filters}.bed",
        bim="{prefix}/{name}{filters}.bim",
        fam="{prefix}/{name}{filters}.fam",
    params:
        out_prefix="{prefix}/{name}{filters}",
    output:
        het="{prefix}/{name}{filters}.het",
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2")
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--het "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"
