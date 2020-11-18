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


################################################################################
# LD Pruning
################################################################################
rule approx_ld:
    """Filter out highly correlated variants within a sliding window.

    Variants in linkage disequilibrium (LD) cannot be separated in GWAS
    analysis. Here we use a 50 variant sliding window to identify variants
    with high LD (i.e., correlation) and prune variants to leave a single
    representative variant for that window. We also filter out variants below
    a MAF threshold.

    .. note::
        `--indep-pairwise <window size>['kb'] <step size (variant ct)> <r^2 threshold>`

        - <window size> sliding window size in variant counts or kilobase
          used to look at correlation
        - <step size> the number of variant counts to slide the window at
          each iteration
        -  <r^2 threshold> pairwise correlation threshold

    .. warning::
        Is the current window size correct? Or should it include the `kb`
        modifier?
    """
    input:
        bed=rules.maf_filter.output.bed,
        bim=rules.maf_filter.output.bim,
        fam=rules.maf_filter.output.fam,
    params:
        in_prefix=rules.maf_filter.params.out_prefix,
        out_prefix="snp_filters/ld_prune/ldPruneList",
        r2=cfg.config.software_params.ld_prune_r2, # r2 threshold: currently 0.1
    output:
        to_keep="snp_filters/ld_prune/ldPruneList.prune.in", # Markers in approx. linkage equilibrium
        to_remove="snp_filters/ld_prune/ldPruneList.prune.out", # Markers in LD
        nosex="snp_filters/ld_prune/ldPruneList.nosex", # Markers in LD
    log:
        "snp_filters/ld_prune/ldPruneList.log",
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
        "--indep-pairwise 50 5 {params.r2}  "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"


rule extract_ld_prune:
    """Subsets the dataset only keeping variants in linkage equilibrium.

    This step keeps only makers in the ``*.prune.in`` and creates a new
    ``plink`` dataset.
    """
    input:
        bed=rules.maf_filter.output.bed,
        bim=rules.maf_filter.output.bim,
        fam=rules.maf_filter.output.fam,
        prune=rules.approx_ld.output.to_keep,
    params:
        in_prefix=rules.maf_filter.params.out_prefix,
        out_prefix="snp_filters/ld_prune/samples",
    output:
        bed="snp_filters/ld_prune/samples.bed",
        bim="snp_filters/ld_prune/samples.bim",
        fam="snp_filters/ld_prune/samples.fam",
        nosex="snp_filters/ld_prune/samples.nosex",
    log:
        "snp_filters/ld_prune/samples.log",
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
        "--extract {input.prune} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"
