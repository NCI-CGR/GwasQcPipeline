"""Filter SNPs and Samples that contain a large number of missing genotype calls.

Call rate filters remove SNPs and samples with a large amount of missingness.
We currently support 2 separate call rate levels. The first level (cr1) is
usually 80% and is used only for reporting. The second level (cr2) is usually
95% and is used to generate a cleaned data set for downstream processing.

The order of filtering (i.e., SNP then Sample or Sample then SNP) will cause
small differences in output. We have opted to minimize the number of Samples
filtered by first filtering SNPs.

.. decision::
    Maximize the number of samples removed due to missingness by first remove
    SNPs with high amounts of missingness.

Input:
    - plink_start/samples.bed
    - plink_start/samples.bim
    - plink_start/samples.fam

Output:

    - Call Rate 1

        - plink_filter_call_rate_1/snps.bed
        - plink_filter_call_rate_1/snps.bim
        - plink_filter_call_rate_1/snps.fam

        - plink_filter_call_rate_1/samples.bed
        - plink_filter_call_rate_1/samples.bim
        - plink_filter_call_rate_1/samples.fam

    - Call Rate 2

        - plink_filter_call_rate_2/snps.bed
        - plink_filter_call_rate_2/snps.bim
        - plink_filter_call_rate_2/snps.fam

        - plink_filter_call_rate_2/samples.bed
        - plink_filter_call_rate_2/samples.bim
        - plink_filter_call_rate_2/samples.fam
"""


################################################################################
# Call Rate 1 (config.software_params.{snp,samp}_cr_1) [default: 0.8]
################################################################################
rule snp_call_rate_filter_1:
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
    group:
        "call_rate_filters"
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


rule sample_call_rate_filter_1:
    input:
        bed=rules.snp_call_rate_filter_1.output.bed,
        bim=rules.snp_call_rate_filter_1.output.bim,
        fam=rules.snp_call_rate_filter_1.output.fam,
    params:
        mind=1 - cfg.config.software_params.samp_cr_1,
        in_prefix=rules.snp_call_rate_filter_1.params.out_prefix,
        out_prefix="plink_filter_call_rate_1/samples",
    output:
        bed="plink_filter_call_rate_1/samples.bed",
        bim="plink_filter_call_rate_1/samples.bim",
        fam="plink_filter_call_rate_1/samples.fam",
    group:
        "call_rate_filters"
    log:
        "plink_filter_call_rate_1/samples.log",
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
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"


################################################################################
# Call Rate 2 (config.software_params.{snp,samp}_cr_2) [default: 0.95]
################################################################################
rule snp_call_rate_filter_2:
    input:
        bed=rules.sample_call_rate_filter_1.output.bed,
        bim=rules.sample_call_rate_filter_1.output.bim,
        fam=rules.sample_call_rate_filter_1.output.fam,
    params:
        geno=1 - cfg.config.software_params.snp_cr_2,
        in_prefix=rules.sample_call_rate_filter_1.params.out_prefix,
        out_prefix="plink_filter_call_rate_2/snps",
    output:
        bed="plink_filter_call_rate_2/snps.bed",
        bim="plink_filter_call_rate_2/snps.bim",
        fam="plink_filter_call_rate_2/snps.fam",
    group:
        "call_rate_filters"
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


rule sample_call_rate_filter_2:
    input:
        bed=rules.snp_call_rate_filter_2.output.bed,
        bim=rules.snp_call_rate_filter_2.output.bim,
        fam=rules.snp_call_rate_filter_2.output.fam,
    params:
        mind=1 - cfg.config.software_params.samp_cr_2,
        in_prefix=rules.snp_call_rate_filter_2.params.out_prefix,
        out_prefix="plink_filter_call_rate_2/samples",
    output:
        bed="plink_filter_call_rate_2/samples.bed",
        bim="plink_filter_call_rate_2/samples.bim",
        fam="plink_filter_call_rate_2/samples.fam",
    group:
        "call_rate_filters"
    log:
        "plink_filter_call_rate_2/samples.log",
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
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"


rule plink_call_rate_stats:
    """Runs ``plink`` missingness statistics.

    .. rubric:: Statistics
        :imiss`: sample-based missing data report
        :lmiss`: variant-based missing data report
    """
    input:
        bed="{prefix}/samples.bed",
        bim="{prefix}/samples.bim",
        fam="{prefix}/samples.fam",
    params:
        in_prefix="{prefix}/samples",
        out_prefix="{prefix}/samples",
    output:
        imiss="{prefix}/samples.imiss",
        lmiss="{prefix}/samples.lmiss",
    group:
        "call_rate_filters"
    resources:
        mem=10000,
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    shell:
        "plink "
        "--bfile {params.in_prefix} "
        "--missing "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"
