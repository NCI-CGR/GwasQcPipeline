from cgr_gwas_qc import load_config

cfg = load_config()


################################################################################
# All Targets
################################################################################
rule all_sample_qc:
    input:
        "sample_level/sample_qc.csv",
        "sample_level/snp_qc.csv",
        "sample_level/concordance/KnownReplicates.csv",
        "sample_level/concordance/InternalQcKnown.csv",
        "sample_level/concordance/StudySampleKnown.csv",
        "sample_level/concordance/UnknownReplicates.csv",
        "sample_level/summary_stats.txt",
        "sample_level/qc_failures/low_call_rate.txt",
        "sample_level/qc_failures/contaminated.txt",
        "sample_level/qc_failures/sex_discordant.txt",
        "sample_level/qc_failures/replicate_discordant.txt",
        "sample_level/internal_controls.txt",
        "sample_level/call_rate.png",
        "sample_level/chrx_inbreeding.png",
        "sample_level/ancestry.png",


################################################################################
# Imports
################################################################################
include: cfg.modules("common")
include: cfg.modules("graf")
include: cfg.modules("plink_filters")
include: cfg.modules("plink_stats")
include: cfg.subworkflow("entry_points")
include: cfg.subworkflow("contamination")


################################################################################
# Workflow Rules
################################################################################
rule update_snps_to_1kg_rsID:
    """Update SNP IDs to rsID from the 1KG project.

    Update study marker IDs to correspond with GRAF's fingerprints which are
    based on the 1KG project rsIDs.
    """
    input:
        bed="{prefix}.bed",
        bim="{prefix}.bim",
        fam="{prefix}.fam",
        vcf=cfg.config.reference_files.thousand_genome_vcf,
    output:
        bed="{prefix}_1kg_rsID.bed",
        bim="{prefix}_1kg_rsID.bim",
        fam="{prefix}_1kg_rsID.fam",
        id_map="{prefix}_1kg_rsID.csv",
    script:
        "../scripts/update_snps_to_1kg_rsID.py"


# -------------------------------------------------------------------------------
# Replicate Concordance
# -------------------------------------------------------------------------------
rule sample_concordance_plink:
    """Summarize sample concordance.

    Splits the concordance table into known and unknown replicates.
    """
    input:
        "sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.concordance.csv".format(
            maf=cfg.config.software_params.maf_for_ibd, ld=cfg.config.software_params.ld_prune_r2,
        ),
    output:
        "sample_level/concordance/plink.csv",
    shell:
        "cp {input[0]} {output[0]}"


rule sample_concordance_graf:
    """Estimate relatedness among samples.

    Outputs a table with pairwise samples and their genotypic relationship.
    GRAF uses a set of 10K pre-screened SNPs, so we can directly use Call
    Rate 2 filtered samples.
    """
    input:
        fpg="sample_level/call_rate_2/samples_1kg_rsID.fpg",
    output:
        "sample_level/concordance/graf.tsv",
    conda:
        cfg.conda("graf")
    log:
        "sample_level/concordance/graf.log",
    shell:
        "graf "
        "-geno {input.fpg} "
        "-type 4 "
        "-out {output[0]} "
        "> {log} 2>&1 "
        "|| exit_code=$?; if [ $exit_code -ne 1 ]; then exit $exit_code; fi" # GRAF returns an exit code of 1, this captures it so snakemake will actually run.


rule sample_concordance_king:
    input:
        bed="sample_level/call_rate_2/samples.bed",
        bim="sample_level/call_rate_2/samples.bim",
        fam="sample_level/call_rate_2/samples.fam",
    output:
        within_family="sample_level/concordance/king.kin",
        between_family="sample_level/concordance/king.kin0",
        within_family_X="sample_level/concordance/kingX.kin",
        between_family_X="sample_level/concordance/kingX.kin0",
        segments="sample_level/concordance/kingallsegs.txt",
    params:
        out_prefix="sample_level/concordance/king",
    conda:
        cfg.conda("king")
    log:
        "sample_level/concordance/king.log",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * 8 * attempt,
        time_hr=lambda wildcards, attempt: 2 * attempt,
    shell:
        "king -b {input.bed} --related --degree 3 --prefix {params.out_prefix} --cpu {threads} > {log} 2>&1 "
        "&& touch {output.within_family} "
        "&& touch {output.between_family} "
        "&& touch {output.within_family_X} "
        "&& touch {output.between_family_X} "
        "&& touch {output.segments} " # king does not always output all files so touch for snakemake


rule sample_concordance:
    input:
        sample_sheet_csv="cgr_sample_sheet.csv",
        plink_file=rules.sample_concordance_plink.output[0],
        graf_file=rules.sample_concordance_graf.output[0],
        king_file=rules.sample_concordance_king.output.between_family,
    output:
        "sample_level/concordance/summary.csv",
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * 12 * attempt,
        time_hr=lambda wildcards, attempt: 2 * attempt,
    script:
        "../scripts/sample_concordance.py"


rule split_sample_concordance:
    input:
        rules.sample_concordance.output[0],
    output:
        known_csv="sample_level/concordance/KnownReplicates.csv",
        known_qc_csv="sample_level/concordance/InternalQcKnown.csv",
        known_study_csv="sample_level/concordance/StudySampleKnown.csv",
        unknown_csv="sample_level/concordance/UnknownReplicates.csv",
    script:
        "../scripts/split_sample_concordance.py"


# -------------------------------------------------------------------------------
# Ancestry
# -------------------------------------------------------------------------------
rule graf_ancestry:
    """Estimate ancestry for each sample."""
    input:
        fpg="sample_level/call_rate_2/samples_1kg_rsID.fpg",
    output:
        "sample_level/ancestry/graf_populations.txt",
    log:
        "sample_level/ancestry/graf_populations.log",
    envmodules:
        cfg.envmodules("graf"),
    conda:
        cfg.conda("graf")
    shell:
        "graf "
        "-geno {input.fpg} "
        "-pop {output[0]} "
        "> {log} 2>&1 "
        "|| exit_code=$?; if [ $exit_code -ne 1 ]; then exit $exit_code; fi" # GRAF returns an exit code of 1, this captures it so snakemake will actually run.


rule graf_ancestry_table:
    """Create summary table with ancestry calls."""
    input:
        rules.graf_ancestry.output[0],
    output:
        "sample_level/ancestry/graf_ancestry.txt",
    envmodules:
        cfg.envmodules("graf"),
    conda:
        cfg.conda("graf")
    shell:
        "PlotPopulations.pl {input[0]} {output[0]} "


# -------------------------------------------------------------------------------
# SNP Summary Table
# -------------------------------------------------------------------------------
rule snp_qc_table:
    input:
        initial="sample_level/samples.lmiss",
        cr1="sample_level/call_rate_1/samples.lmiss",
        cr2="sample_level/call_rate_2/samples.lmiss",
        thousand_genomes="sample_level/call_rate_2/samples_1kg_rsID.csv",
    output:
        "sample_level/snp_qc.csv",
    script:
        "../scripts/snp_qc_table.py"


# -------------------------------------------------------------------------------
# Sample Summary Table
# -------------------------------------------------------------------------------
def _contam(wildcards):
    uf, wf = cfg.config.user_files, cfg.config.workflow_params
    if uf.idat_pattern and uf.gtc_pattern and wf.remove_contam:
        return "sample_level/contamination.csv"
    return []


def _intensity(wildcards):
    uf, wf = cfg.config.user_files, cfg.config.workflow_params
    if uf.idat_pattern and uf.gtc_pattern and wf.remove_contam:
        return "sample_level/median_idat_intensity.csv"
    return []


rule sample_qc_table:
    input:
        sample_sheet_csv="cgr_sample_sheet.csv",
        imiss_start="sample_level/samples.imiss",
        imiss_cr1="sample_level/call_rate_1/samples.imiss",
        imiss_cr2="sample_level/call_rate_2/samples.imiss",
        sexcheck_cr1="sample_level/call_rate_1/samples.sexcheck",
        ancestry="sample_level/ancestry/graf_ancestry.txt",
        sample_concordance_csv="sample_level/concordance/summary.csv",
        contam=_contam,
        intensity=_intensity,
    params:
        remove_contam=cfg.config.workflow_params.remove_contam,
        remove_rep_discordant=cfg.config.workflow_params.remove_rep_discordant,
    output:
        "sample_level/sample_qc.csv",
    script:
        "../scripts/sample_qc_table.py"


rule sample_qc_summary_stats:
    input:
        rules.sample_qc_table.output[0],
    output:
        "sample_level/summary_stats.txt",
    script:
        "../scripts/sample_qc_summary_stats.py"


rule sample_lists_from_qc_flags:
    input:
        all_samples=rules.sample_qc_table.output[0],
    output:
        cr="sample_level/qc_failures/low_call_rate.txt",
        contam="sample_level/qc_failures/contaminated.txt",
        sex="sample_level/qc_failures/sex_discordant.txt",
        rep="sample_level/qc_failures/replicate_discordant.txt",
        ctrl="sample_level/internal_controls.txt",
    script:
        "../scripts/sample_lists_from_qc_flags.py"


# -------------------------------------------------------------------------------
# Plots
# -------------------------------------------------------------------------------
rule plot_call_rate:
    input:
        sample_qc=rules.sample_qc_table.output[0],
        snp_qc=rules.snp_qc_table.output[0],
    params:
        sample_cr1=cfg.config.software_params.sample_call_rate_1,
        snp_cr1=cfg.config.software_params.snp_call_rate_1,
        sample_cr2=cfg.config.software_params.sample_call_rate_2,
        snp_cr2=cfg.config.software_params.snp_call_rate_2,
    output:
        "sample_level/call_rate.png",
    script:
        "../scripts/plot_call_rate.py"


rule plot_chrx_inbreeding:
    input:
        rules.sample_qc_table.output[0],
    output:
        "sample_level/chrx_inbreeding.png",
    script:
        "../scripts/plot_chrx_inbreeding.py"


rule plot_ancestry:
    input:
        rules.sample_qc_table.output[0],
    output:
        "sample_level/ancestry.png",
    script:
        "../scripts/plot_ancestry.py"
