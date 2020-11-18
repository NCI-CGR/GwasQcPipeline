"""Module containing SNP filters"""


################################################################################
# MAF Filtering
################################################################################
rule maf_filter:
    """Filter SNPs based on minor allele frequency.

    .. warning::
        The current default MAF threshold is 0.2. This seems really high.
    """
    input:
        bed="plink_filter_call_rate_2/samples.bed",
        bim="plink_filter_call_rate_2/samples.bim",
        fam="plink_filter_call_rate_2/samples.fam",
    params:
        in_prefix="plink_filter_call_rate_2/samples",
        out_prefix="snp_filters/maf/samples",
        maf=cfg.config.software_params.maf_for_ibd,
    output:
        bed="snp_filters/maf/samples.bed",
        bim="snp_filters/maf/samples.bim",
        fam="snp_filters/maf/samples.fam",
        nosex="snp_filters/maf/samples.nosex",
    log:
        "snp_filters/maf/samples.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 2
    resources:
        mem=10000,
    shell:
        "plink "
        "--bfile {params.in_prefix} "
        "--maf {params.maf} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"
