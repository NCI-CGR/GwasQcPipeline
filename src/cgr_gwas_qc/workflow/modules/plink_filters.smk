include: cfg.modules("common.smk")


rule sample_call_rate_1:
    input:
        bed="sample_level/samples.bed",
        bim="sample_level/samples.bim",
        fam="sample_level/samples.fam",
    params:
        mind=1 - cfg.config.software_params.sample_call_rate_1,
        out_prefix="sample_level/call_rate_1/samples_p1",
    output:
        bed=temp("sample_level/call_rate_1/samples_p1.bed"),
        bim=temp("sample_level/call_rate_1/samples_p1.bim"),
        fam=temp("sample_level/call_rate_1/samples_p1.fam"),
        nosex=temp("sample_level/call_rate_1/samples_p1.nosex"),
    log:
        "sample_level/call_rate_1/samples_p1.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
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


rule snp_call_rate_1:
    input:
        bed=rules.sample_call_rate_1.output.bed,
        bim=rules.sample_call_rate_1.output.bim,
        fam=rules.sample_call_rate_1.output.fam,
    params:
        geno=1 - cfg.config.software_params.snp_call_rate_1,
        out_prefix="sample_level/call_rate_1/samples",
    output:
        bed="sample_level/call_rate_1/samples.bed",
        bim="sample_level/call_rate_1/samples.bim",
        fam="sample_level/call_rate_1/samples.fam",
        nosex="sample_level/call_rate_1/samples.nosex",
    log:
        "sample_level/call_rate_1/samples.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
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


rule sample_call_rate_2:
    input:
        bed="sample_level/call_rate_1/samples.bed",
        bim="sample_level/call_rate_1/samples.bim",
        fam="sample_level/call_rate_1/samples.fam",
    params:
        mind=1 - cfg.config.software_params.sample_call_rate_2,
        out_prefix="sample_level/call_rate_2/samples_p1",
    output:
        bed=temp("sample_level/call_rate_2/samples_p1.bed"),
        bim=temp("sample_level/call_rate_2/samples_p1.bim"),
        fam=temp("sample_level/call_rate_2/samples_p1.fam"),
        nosex=temp("sample_level/call_rate_2/samples_p1.nosex"),
    log:
        "sample_level/call_rate_2/samples_p1.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
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


rule snp_call_rate_2:
    input:
        bed=rules.sample_call_rate_2.output.bed,
        bim=rules.sample_call_rate_2.output.bim,
        fam=rules.sample_call_rate_2.output.fam,
    params:
        geno=1 - cfg.config.software_params.snp_call_rate_2,
        out_prefix="sample_level/call_rate_2/samples",
    output:
        bed="sample_level/call_rate_2/samples.bed",
        bim="sample_level/call_rate_2/samples.bim",
        fam="sample_level/call_rate_2/samples.fam",
        nosex="sample_level/call_rate_2/samples.nosex",
    log:
        "sample_level/call_rate_2/samples.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
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
    """Filter SNPs based on minor allele frequency.

    .. warning::
        The current default MAF threshold is 0.2. This seems really high.
    """
    input:
        bed="{prefix}/{name}{filters}.bed",
        bim="{prefix}/{name}{filters}.bim",
        fam="{prefix}/{name}{filters}.fam",
    params:
        maf="{maf}",
        out_prefix="{prefix}/{name}{filters}_maf{maf}",
    output:
        bed=temp("{prefix}/{name}{filters}_maf{maf}.bed"),
        bim=temp("{prefix}/{name}{filters}_maf{maf}.bim"),
        fam=temp("{prefix}/{name}{filters}_maf{maf}.fam"),
        nosex=temp("{prefix}/{name}{filters}_maf{maf}.nosex"),
    log:
        "{prefix}/{name}{filters}_maf{maf}.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
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
        bed="{prefix}/{name}{filters}.bed",
        bim="{prefix}/{name}{filters}.bim",
        fam="{prefix}/{name}{filters}.fam",
    params:
        r2="{ld}", # r2 threshold: currently 0.1
        out_prefix="{prefix}/{name}{filters}_ld{ld}",
    output:
        to_keep=temp("{prefix}/{name}{filters}_ld{ld}.prune.in"), # Markers in approx. linkage equilibrium
        to_remove=temp("{prefix}/{name}{filters}_ld{ld}.prune.out"), # Markers in LD
        nosex=temp("{prefix}/{name}{filters}_ld{ld}.nosex"), # Markers in LD
    log:
        "{prefix}/{name}{filters}_ld{ld}.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
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


rule ld_prune:
    """Subsets the dataset only keeping variants in linkage equilibrium.

    This step keeps only makers in the ``*.prune.in`` and creates a new
    ``plink`` dataset.
    """
    input:
        bed="{prefix}/{name}{filters}.bed",
        bim="{prefix}/{name}{filters}.bim",
        fam="{prefix}/{name}{filters}.fam",
        prune="{prefix}/{name}{filters}_ld{ld}.prune.in",
    params:
        out_prefix="{prefix}/{name}{filters}_ld{ld}_pruned",
    output:
        bed="{prefix}/{name}{filters}_ld{ld}_pruned.bed",
        bim="{prefix}/{name}{filters}_ld{ld}_pruned.bim",
        fam="{prefix}/{name}{filters}_ld{ld}_pruned.fam",
        nosex="{prefix}/{name}{filters}_ld{ld}_pruned.nosex",
    log:
        "{prefix}/{name}{filters}_ld{ld}_pruned.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--extract {input.prune} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule snps_only_filter:
    """Exclude all variants with one or more multi-character allele codes"""
    input:
        bed="{prefix}/{name}{filters}.bed",
        bim="{prefix}/{name}{filters}.bim",
        fam="{prefix}/{name}{filters}.fam",
    params:
        out_prefix="{prefix}/{name}{filters}_snps",
    output:
        bed=temp("{prefix}/{name}{filters}_snps.bed"),
        bim=temp("{prefix}/{name}{filters}_snps.bim"),
        fam=temp("{prefix}/{name}{filters}_snps.fam"),
        nosex=temp("{prefix}/{name}{filters}_snps.nosex"),
    log:
        "{prefix}/{name}{filters}_snps.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
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
        bed="{prefix}/{name}{filters}.bed",
        bim="{prefix}/{name}{filters}.bim",
        fam="{prefix}/{name}{filters}.fam",
    params:
        out_prefix="{prefix}/{name}{filters}_autosome",
    output:
        bed=temp("{prefix}/{name}{filters}_autosome.bed"),
        bim=temp("{prefix}/{name}{filters}_autosome.bim"),
        fam=temp("{prefix}/{name}{filters}_autosome.fam"),
        nosex=temp("{prefix}/{name}{filters}_autosome.nosex"),
    log:
        "{prefix}/{name}{filters}_autosome.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
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


rule cleaned_filter:
    """Tell snakemake to keep the file.

    The link filter rules are designed to be strung together. Intermediate
    files are automatically deleted by snakemake. This "filter" is used to
    create files that are kept by snakemake.
    """
    input:
        bed="{prefix}/{name}{filters}.bed",
        bim="{prefix}/{name}{filters}.bim",
        fam="{prefix}/{name}{filters}.fam",
    params:
        out_prefix="{prefix}/{name}{filters}_cleaned",
    output:
        bed="{prefix}/{name}{filters}_cleaned.bed",
        bim="{prefix}/{name}{filters}_cleaned.bim",
        fam="{prefix}/{name}{filters}_cleaned.fam",
        nosex="{prefix}/{name}{filters}_cleaned.nosex",
    log:
        "{prefix}/{name}{filters}_cleaned.log",
    envmodules:
        cfg.envmodules("plink2"),
    conda:
        cfg.conda("plink2.yml")
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
