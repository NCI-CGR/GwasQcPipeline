def _snp_call_rate_filter(wildcards):
    if wildcards.cr == "2":
        return {
            "bed": "{prefix}/call_rate_1/samples.bed",
            "bim": "{prefix}/call_rate_1/samples.bim",
            "fam": "{prefix}/call_rate_1/samples.fam",
        }

    return {
        "bed": "{prefix}/samples.bed",
        "bim": "{prefix}/samples.bim",
        "fam": "{prefix}/samples.fam",
    }


rule snp_call_rate_filter:
    input:
        unpack(_snp_call_rate_filter),
    params:
        geno=lambda wc: 1 - float(cfg.config.software_params.dict().get(f"snp_cr_{wc.cr}")),
        out_prefix="{prefix}/call_rate_{cr}/snps",
    output:
        bed=temp("{prefix}/call_rate_{cr}/snps.bed"),
        bim=temp("{prefix}/call_rate_{cr}/snps.bim"),
        fam=temp("{prefix}/call_rate_{cr}/snps.fam"),
    wildcard_constraints:
        snp_cr="[1, 2]",
    log:
        "{prefix}/call_rate_{cr}/snps.log",
    group:
        "call_rate_filters"
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
        "--geno {params.geno} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"


rule sample_call_rate_filter:
    input:
        bed="{prefix}/call_rate_{cr}/snps.bed",
        bim="{prefix}/call_rate_{cr}/snps.bim",
        fam="{prefix}/call_rate_{cr}/snps.fam",
    params:
        mind=lambda wc: 1 - float(cfg.config.software_params.dict().get(f"samp_cr_{wc.cr}")),
        out_prefix="{prefix}/call_rate_{cr}/samples",
    output:
        bed=temp("{prefix}/call_rate_{cr}/samples.bed"),
        bim=temp("{prefix}/call_rate_{cr}/samples.bim"),
        fam=temp("{prefix}/call_rate_{cr}/samples.fam"),
    wildcard_constraints:
        samp_cr="[1, 2]",
    log:
        "{prefix}/call_rate_{cr}/samples.log",
    group:
        "call_rate_filters"
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
        "--mind {params.mind} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"


rule maf_filter:
    """Filter SNPs based on minor allele frequency.

    .. warning::
        The current default MAF threshold is 0.2. This seems really high.
    """
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        maf=lambda wc: float(wc.maf),
        out_prefix="{prefix}_maf{maf}",
    output:
        bed=temp("{prefix}_maf{maf}.bed"),
        bim=temp("{prefix}_maf{maf}.bim"),
        fam=temp("{prefix}_maf{maf}.fam"),
        nosex=temp("{prefix}_maf{maf}.nosex"),
    wildcard_constraints:
        maf="0\.\d+",
    log:
        "{prefix}_maf{maf}.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 2
    resources:
        mem=10000,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--maf {params.maf} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"


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
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
    params:
        r2=lambda wc: float(wc.ld), # r2 threshold: currently 0.1
        out_prefix="{prefix}_ld{ld}",
    output:
        to_keep=temp("{prefix}_ld{ld}.prune.in"), # Markers in approx. linkage equilibrium
        to_remove=temp("{prefix}_ld{ld}.prune.out"), # Markers in LD
        nosex=temp("{prefix}_ld{ld}.nosex"), # Markers in LD
    wildcard_constraints:
        ld="0\.\d+",
    log:
        "{prefix}_ld{ld}.log",
    group:
        "ld"
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 2
    resources:
        mem=10000,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--indep-pairwise 50 5 {params.r2}  "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"


rule ld_prune:
    """Subsets the dataset only keeping variants in linkage equilibrium.

    This step keeps only makers in the ``*.prune.in`` and creates a new
    ``plink`` dataset.
    """
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
        prune="{prefix}_ld{ld}.prune.in",
    params:
        out_prefix="{prefix}_ld{ld}_pruned",
    output:
        bed="{prefix}_ld{ld}_pruned.bed",
        bim="{prefix}_ld{ld}_pruned.bim",
        fam="{prefix}_ld{ld}_pruned.fam",
        nosex="{prefix}_ld{ld}_pruned.nosex",
    wildcard_constraints:
        ld="0\.\d+",
    log:
        "{prefix}_ld{ld}_pruned.log",
    group:
        "ld"
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 2
    resources:
        mem=10000,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--extract {input.prune} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
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
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 2
    resources:
        mem=10000,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--snps-only "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
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
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 2
    resources:
        mem=10000,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--autosome "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"


rule cleaned_filter:
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
        out_prefix="{prefix}_cleaned",
    output:
        bed="{prefix}_cleaned.bed",
        bim="{prefix}_cleaned.bim",
        fam="{prefix}_cleaned.fam",
        nosex="{prefix}_cleaned.nosex",
    log:
        "{prefix}_cleaned.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: 2
    resources:
        mem=10000,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix}"
