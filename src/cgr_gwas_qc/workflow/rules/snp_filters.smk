"""Marker level filters"""

################################################################################
# SNP Call Rate Filters (missingness)
################################################################################
"""Removes SNPs with low genotype call rates.

This module uses ``plink`` to filter SNPs that have low genotype
call rates. This is an indication of problems ranging from low DNA quality or
possibly a sample from a different population. We use two successive call
rates from the ``config.yml``.

- ``software_params.snp_cr_1``
- ``software_params.snp_cr_2``

Outputs are used in Sample level call rate filters.

Inputs:

- Call Rate 1
    - plink_start/samples.bed
    - plink_start/samples.bim
    - plink_start/samples.fam

- Call Rate 2
    - plink_filter_call_rate_2/snps.bed
    - plink_filter_call_rate_2/snps.bim
    - plink_filter_call_rate_2/snps.fam

Outputs:

- Call Rate 1
    - plink_filter_call_rate_1/snps.bed
    - plink_filter_call_rate_1/snps.bim
    - plink_filter_call_rate_1/snps.fam

- Call Rate 2
    - plink_filter_call_rate_2/snps.bed
    - plink_filter_call_rate_2/snps.bim
    - plink_filter_call_rate_2/snps.fam
"""


rule snp_call_rate_filter_1:
    """Removes SNPs below call rate filter 1 and outputs a binary ``plink`` format.

    Call rate filter 1 is usually (1 - 0.8) missing for samples and variants.
    """
    input:
        bed="plink_start/samples.bed",
        bim="plink_start/samples.bim",
        fam="plink_start/samples.fam",
    params:
        geno=1 - cfg.config.software_params.snp_cr_1,
        in_prefix="plink_start/samples",
        out_prefix="plink_filter_call_rate_1/snps",
    output:
        bed="plink_filter_call_rate_1/snps.bed",
        bim="plink_filter_call_rate_1/snps.bim",
        fam="plink_filter_call_rate_1/snps.fam",
    log:
        "plink_filter_call_rate_1/snps.log",
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
        "--geno {params.geno} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"


rule snp_call_rate_filter_2:
    """Removes SNPs below call rate filter 1 and outputs a binary ``plink`` format.

    Call rate filter 2 is usually (1 - 0.95) missing for samples and variants.
    """
    input:
        bed="plink_filter_call_rate_1/samples.bed",
        bim="plink_filter_call_rate_1/samples.bim",
        fam="plink_filter_call_rate_1/samples.fam",
    params:
        geno=1 - cfg.config.software_params.snp_cr_2,
        in_prefix="plink_filter_call_rate_1/samples",
        out_prefix="plink_filter_call_rate_2/snps",
    output:
        bed="plink_filter_call_rate_2/snps.bed",
        bim="plink_filter_call_rate_2/snps.bim",
        fam="plink_filter_call_rate_2/snps.fam",
    log:
        "plink_filter_call_rate_2/snps.log",
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
        "--geno {params.geno} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"
