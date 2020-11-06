"""Removes Samples and SNPs with low genotype call rates.

This module uses ``plink`` to filter Samples and SNPs that have low genotype
call rates. This is an indication of problems ranging from low DNA quality or
possibly a sample from a different population. We use two successive call
rates from the ``config.yml``.

- ``software_params.samp_cr_1``, ``software_params.snp_cr_1``
- ``software_params.samp_cr_2``, ``software_params.snp_cr_2``

Downstream QC steps use the Call Rate 2 outputs, but Call Rate 1 results are
included in final reports.

Inputs:

- plink_start/samples.bed
- plink_start/samples.bim
- plink_start/samples.fam

Outputs:

- Call Rate 1
    - plink_filter_call_rate_1/samples.bed
    - plink_filter_call_rate_1/samples.bim
    - plink_filter_call_rate_1/samples.fam

- Call Rate 2
    - plink_filter_call_rate_2/samples.bed
    - plink_filter_call_rate_2/samples.bim
    - plink_filter_call_rate_2/samples.fam
"""
rule call_rate_filter_1:
    """Removes samples and SNPs below call rate filter 1 and outputs a binary ``plink`` format.

    Call rate filter 1 is usually (1 - 0.8) missing for samples and variants.
    """
    input:
        bed="plink_start/samples.bed",
        bim="plink_start/samples.bim",
        fam="plink_start/samples.fam",
    params:
        mind=1 - cfg.config.software_params.samp_cr_1,
        geno=1 - cfg.config.software_params.snp_cr_1,
        in_prefix="plink_start/samples",
        out_prefix="plink_filter_call_rate_1/samples",
    output:
        bed="plink_filter_call_rate_1/samples.bed",
        bim="plink_filter_call_rate_1/samples.bim",
        fam="plink_filter_call_rate_1/samples.fam",
    log: "plink_filter_call_rate_1/samples.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 20
    resources:
        mem=10000,
    shell:
        "plink "
        "--bfile {params.in_prefix} "
        "--mind {params.mind} "
        "--geno {params.geno} "
        "--make-bed "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"


rule call_rate_filter_2:
    """Removes samples and SNPs below call rate filter 2 and outputs a binary ``plink`` format.

    Call rate filter 2 is usually (1 - 0.95) missing for samples and variants.
    """
    input:
        bed="plink_filter_call_rate_1/samples.bed",
        bim="plink_filter_call_rate_1/samples.bim",
        fam="plink_filter_call_rate_1/samples.fam",
    params:
        mind=1 - cfg.config.software_params.samp_cr_2,
        geno=1 - cfg.config.software_params.snp_cr_2,
        in_prefix="plink_filter_call_rate_1/samples",
        out_prefix="plink_filter_call_rate_2/samples",
    output:
        bed="plink_filter_call_rate_2/samples.bed",
        bim="plink_filter_call_rate_2/samples.bim",
        fam="plink_filter_call_rate_2/samples.fam",
    log: "plink_filter_call_rate_2/samples.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 20
    resources:
        mem=10000,
    shell:
        "plink "
        "--bfile {params.in_prefix} "
        "--mind {params.mind} "
        "--geno {params.geno} "
        "--make-bed "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"
