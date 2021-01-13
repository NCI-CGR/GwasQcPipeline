"""Module containing SNP filters"""


################################################################################
# Filter SNPs not in 1000 Genome (1KG)
################################################################################
rule thousG_match:
    """Flags SNPs and flips alleles to match the 1KG project.

    This is a custom script that compares markers with the 1KG vcf file. For
    a given SNP, if 1KG uses the compliment this script will flip the allele.
    It then flags SNPs for removal if they:

    1. Don't match a SNP in 1KG
    2. Ambiguous alleles (i.e., A/T or C/G)
    3. Indels
    4. SNPs on chromosomes other than 1-22 and X (converts 23 to X)
    """
    input:
        bim=rules.extract_ld_prune.output.bim,
        vcf=cfg.config.reference_files.thousand_genome_vcf,
    output:
        snps_to_remove="snp_filters/plink_thousG_match/thousG_snps_to_remove.txt",
        bim="snp_filters/plink_thousG_match/samples.vcfStrand.bim",
    log:
        "snp_filters/plink_thousG_match/thousG_match.log",
    script:
        "../scripts/bim_filter_vcf.py"


rule remove_not_thousG:
    """Removes SNPs not matching the 1KG project.

    Outputs a new set of ``plink`` files BED + BIM + FAM removing the SNPs
    (and Indels) flagged as not in the 1KG project.
    """
    input:
        bim=rules.thousG_match.output.bim,
        bed=rules.extract_ld_prune.output.bed,
        fam=rules.extract_ld_prune.output.fam,
        snps_to_remove=rules.thousG_match.output.snps_to_remove,
    params:
        out_prefix="snp_filters/plink_thousG_match/samples",
    output:
        "snp_filters/plink_thousG_match/samples.bed",
        "snp_filters/plink_thousG_match/samples.bim",
        "snp_filters/plink_thousG_match/samples.fam",
    log:
        "snp_filters/plink_thousG_match/samples.log",
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
        "--exclude {input.snps_to_remove} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem} "
        "--out {params.out_prefix} "
