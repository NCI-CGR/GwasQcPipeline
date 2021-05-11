from cgr_gwas_qc import load_config

cfg = load_config()


################################################################################
# Filters
################################################################################
rule sample_call_rate_filter:
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        mind=0.95,
        out_prefix="{prefix}.sample_cr_filtered",
    output:
        bed=temp("{prefix}.sample_cr_filtered.bed"),
        bim=temp("{prefix}.sample_cr_filtered.bim"),
        fam=temp("{prefix}.sample_cr_filtered.fam"),
        nosex=temp("{prefix}.sample_cr_filtered.nosex"),
    log:
        "{prefix}.sample_cr_filtered.log",
    conda:
        cfg.conda("plink2")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        tim_hr=lambda wildcards, attempt: attempt ** 2,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--mind {params.mind} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule snp_call_rate_filter:
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        geno=0.95,
        out_prefix="{prefix}.snp_cr_filtered",
    output:
        bed=temp("{prefix}.snp_cr_filtered.bed"),
        bim=temp("{prefix}.snp_cr_filtered.bim"),
        fam=temp("{prefix}.snp_cr_filtered.fam"),
        nosex=temp("{prefix}.snp_cr_filtered.nosex"),
    log:
        "{prefix}.snp_cr_filtered.log",
    conda:
        cfg.conda("plink2")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        tim_hr=lambda wildcards, attempt: attempt ** 2,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--geno {params.geno} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule maf_filter:
    """Filter SNPs based on minor allele frequency."""
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        maf="{maf}",
        out_prefix="{prefix}_maf{maf}",
    output:
        bed=temp("{prefix}_maf{maf}.bed"),
        bim=temp("{prefix}_maf{maf}.bim"),
        fam=temp("{prefix}_maf{maf}.fam"),
        nosex=temp("{prefix}_maf{maf}.nosex"),
    log:
        "{prefix}_maf{maf}.log",
    wildcard_constraints:
        maf="0.\d+",
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
        "--maf {params.maf} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule ld_filter:
    """Subsets the dataset only keeping variants in linkage equilibrium.

    This step keeps only makers in the ``*.prune.in`` and creates a new
    ``plink`` dataset.
    """
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
        to_keep="{prefix}_ld{ld}.prune.in",
    params:
        out_prefix="{prefix}_ld{ld}_pruned",
    output:
        bed=temp("{prefix}_ld{ld}_pruned.bed"),
        bim=temp("{prefix}_ld{ld}_pruned.bim"),
        fam=temp("{prefix}_ld{ld}_pruned.fam"),
        nosex=temp("{prefix}_ld{ld}_pruned.nosex"),
    log:
        "{prefix}_ld{ld}_pruned.log",
    wildcard_constraints:
        ld="0.\d+",
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
        "--extract {input.to_keep} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule snps_only_filter:
    """Exclude all variants with one or more multi-character allele codes"""
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        out_prefix="{prefix}_snps",
    output:
        bed=temp("{prefix}_snps.bed"),
        bim=temp("{prefix}_snps.bim"),
        fam=temp("{prefix}_snps.fam"),
        nosex=temp("{prefix}_snps.nosex"),
    log:
        "{prefix}_snps.log",
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
        "--snps-only "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule autosome_only_filter:
    """Exclude all unplaced and non-autosomal variants"""
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        out_prefix="{prefix}_autosome",
    output:
        bed=temp("{prefix}_autosome.bed"),
        bim=temp("{prefix}_autosome.bim"),
        fam=temp("{prefix}_autosome.fam"),
        nosex=temp("{prefix}_autosome.nosex"),
    log:
        "{prefix}_autosome.log",
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
        "--autosome "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule keep_ids:
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
        to_keep="{prefix}_to_keep.txt",
    params:
        out_prefix="{prefix}.filtered",
    output:
        bed="{prefix}.filtered.bed",
        bim="{prefix}.filtered.bim",
        fam="{prefix}.filtered.fam",
        nosex="{prefix}.filtered.nosex",
    log:
        "{prefix}.filtered.log",
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
        "--keep {input.to_keep} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule remove_ids:
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
        to_remove="{prefix}_to_remove.txt",
    params:
        out_prefix="{prefix}.filtered",
    output:
        bed="{prefix}.filtered.bed",
        bim="{prefix}.filtered.bim",
        fam="{prefix}.filtered.fam",
        nosex="{prefix}.filtered.nosex",
    log:
        "{prefix}.filtered.log",
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
        "--remove {input.to_remove} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule keep_bfile:
    """Tell snakemake to keep the file.

    The link filter rules are designed to be strung together. Intermediate
    files are automatically deleted by snakemake. This "filter" is used to
    create files that are kept by snakemake.
    """
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        out_prefix="{prefix}.keep",
    output:
        bed="{prefix}.keep.bed",
        bim="{prefix}.keep.bim",
        fam="{prefix}.keep.fam",
        nosex="{prefix}.keep.nosex",
    log:
        "{prefix}.keep.log",
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
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


################################################################################
# Converters
################################################################################
rule rename_ids:
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
        id_map="{prefix}.id_map.txt",
    params:
        out_prefix="{prefix}.renamed",
    output:
        bed="{prefix}.renamed.bed",
        bim="{prefix}.renamed.bim",
        fam="{prefix}.renamed.fam",
        nosex="{prefix}.renamed.nosex",
    log:
        "{prefix}.renamed.log",
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
        "--update-ids {input.id_map} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule bed_to_ped:
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    output:
        ped="{prefix}.ped",
        map_="{prefix}.map",
    params:
        out_prefix="{prefix}",
    log:
        "{prefix}.log",
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
        "--out {params.out_prefix}"


################################################################################
# Stats
################################################################################
rule ld:
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
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        r2="{ld}",  # r2 threshold: currently 0.1
        out_prefix="{prefix}_ld{ld}",
    output:
        to_keep=temp("{prefix}_ld{ld}.prune.in"),  # Markers in approx. linkage equilibrium
        to_remove=temp("{prefix}_ld{ld}.prune.out"),  # Markers in LD
        nosex=temp("{prefix}_ld{ld}.nosex"),  # Markers in LD
    log:
        "{prefix}_ld{ld}.log",
    wildcard_constraints:
        ld="0.\d+",
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
        "--indep-pairwise 50 5 {params.r2}  "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule miss:
    """Runs ``plink`` missingness statistics.

    .. rubric:: Statistics
        :imiss`: sample-based missing data report
        :lmiss`: variant-based missing data report
    """
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        out_prefix="{prefix}",
    output:
        imiss="{prefix}.imiss",
        lmiss="{prefix}.lmiss",
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


rule sexcheck:
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        out_prefix="{prefix}",
    output:
        "{prefix}.sexcheck",
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
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


rule frq:
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        out_prefix="{prefix}",
    output:
        "{prefix}.frq",
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
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


rule hwe:
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        out_prefix="{prefix}",
    output:
        "{prefix}.hwe",
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
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


rule genome:
    """Calculates Identity-by-descent.

    This step is trying to find relationships among the samples by
    calculating IBS/IBD. This calculation is no LD-aware so we are doing the
    LD pruning before. This calculation excludes non-autosomes.
    """
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        ibd_min=cfg.config.software_params.ibd_pi_hat_min,
        ibd_max=cfg.config.software_params.ibd_pi_hat_max,
        out_prefix="{prefix}",
    output:
        "{prefix}.genome",
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
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


rule het:
    """Calculates autosomal heterozygosity."""
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        out_prefix="{prefix}",
    output:
        "{prefix}.het",
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
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
