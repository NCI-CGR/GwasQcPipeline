from cgr_gwas_qc import load_config

cfg = load_config()

use_contamination = (
    cfg.config.user_files.idat_pattern
    and cfg.config.user_files.gtc_pattern
    and cfg.config.workflow_params.remove_contam
)


localrules:
    all_sample_qc,
    sample_concordance_plink,
    sample_concordance_summary,
    split_sample_concordance,
    snp_qc_table,
    sample_level_sexcheck,
    sample_qc_table,
    sample_qc_summary_stats,
    sample_lists_from_qc_flags,
    plot_call_rate,
    plot_chrx_inbreeding,
    plot_ancestry,


################################################################################
# Sample QC Targets
################################################################################
targets = [
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
]

if use_contamination:
    targets.append("sample_level/contamination/summary.csv")


rule all_sample_qc:
    input:
        targets,


################################################################################
# Imports
################################################################################
module thousand_genomes:
    snakefile:
        cfg.modules("thousand_genomes")
    config:
        {}


module plink:
    snakefile:
        cfg.modules("plink")


module graf:
    snakefile:
        cfg.modules("graf")


################################################################################
# Workflow Rules
################################################################################
# -------------------------------------------------------------------------------
# Call Rate Filters
# -------------------------------------------------------------------------------
use rule sample_call_rate_filter from plink as sample_call_rate_filter_1 with:
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


use rule snp_call_rate_filter from plink as snp_call_rate_filter_1 with:
    input:
        bed=rules.sample_call_rate_filter_1.output.bed,
        bim=rules.sample_call_rate_filter_1.output.bim,
        fam=rules.sample_call_rate_filter_1.output.fam,
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


use rule sample_call_rate_filter from plink as sample_call_rate_filter_2 with:
    input:
        bed=rules.snp_call_rate_filter_1.output.bed,
        bim=rules.snp_call_rate_filter_1.output.bim,
        fam=rules.snp_call_rate_filter_1.output.fam,
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


use rule snp_call_rate_filter from plink as snp_call_rate_filter_2 with:
    input:
        bed=rules.sample_call_rate_filter_2.output.bed,
        bim=rules.sample_call_rate_filter_2.output.bim,
        fam=rules.sample_call_rate_filter_2.output.fam,
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


# -------------------------------------------------------------------------------
# Call Rate Statistics
# -------------------------------------------------------------------------------
use rule miss from plink as plink_call_rate_initial with:
    input:
        bed="sample_level/samples.bed",
        bim="sample_level/samples.bim",
        fam="sample_level/samples.fam",
    params:
        out_prefix="sample_level/samples",
    output:
        imiss="sample_level/samples.imiss",
        lmiss="sample_level/samples.lmiss",


use rule miss from plink as plink_call_rate_post1 with:
    input:
        bed=rules.snp_call_rate_filter_1.output.bed,
        bim=rules.snp_call_rate_filter_1.output.bim,
        fam=rules.snp_call_rate_filter_1.output.fam,
    params:
        out_prefix="sample_level/call_rate_1/samples",
    output:
        imiss="sample_level/call_rate_1/samples.imiss",
        lmiss="sample_level/call_rate_1/samples.lmiss",


use rule miss from plink as plink_call_rate_post2 with:
    input:
        bed=rules.snp_call_rate_filter_2.output.bed,
        bim=rules.snp_call_rate_filter_2.output.bim,
        fam=rules.snp_call_rate_filter_2.output.fam,
    params:
        out_prefix="sample_level/call_rate_2/samples",
    output:
        imiss="sample_level/call_rate_2/samples.imiss",
        lmiss="sample_level/call_rate_2/samples.lmiss",


# -------------------------------------------------------------------------------
# Contamination
# -------------------------------------------------------------------------------
if use_contamination:

    localrules:
        sample_contamination_verifyIDintensity,

    rule sample_contamination_verifyIDintensity:
        """Aggregate sample contamination scores.

        Aggregates sample level contamination scores into a single file (each
        row is a sample). The script sets ``%Mix`` to ``NA`` if the intensity
        is below the threshold and the file is not in the ``imiss`` file.
        """
        input:
            contamination_file="sample_level/contamination/verifyIDintensity.csv",
            median_intensity_file="sample_level/contamination/median_idat_intensity.csv",
            imiss_file=rules.plink_call_rate_post2.output.imiss,
        params:
            intensity_threshold=cfg.config.software_params.intensity_threshold,
            contam_threshold=cfg.config.software_params.contam_threshold,
        output:
            "sample_level/contamination/summary.csv",
        script:
            "../scripts/agg_contamination.py"


# -------------------------------------------------------------------------------
# Update GSA to 1KG rsID
# -------------------------------------------------------------------------------
use rule update_snps_to_1kg_rsID from thousand_genomes as update_samples_to_1kg_rsIDs with:
    input:
        bed=rules.snp_call_rate_filter_2.output.bed,
        bim=rules.snp_call_rate_filter_2.output.bim,
        fam=rules.snp_call_rate_filter_2.output.fam,
        vcf=cfg.config.reference_files.thousand_genome_vcf,
    output:
        bed="sample_level/call_rate_2/samples_1kg_rsID.bed",
        bim="sample_level/call_rate_2/samples_1kg_rsID.bim",
        fam="sample_level/call_rate_2/samples_1kg_rsID.fam",
        id_map="sample_level/call_rate_2/samples_1kg_rsID.csv",


# -------------------------------------------------------------------------------
# Replicate Concordance
# -------------------------------------------------------------------------------
use rule maf_filter from plink as sample_level_maf_filter with:
    input:
        bed=rules.snp_call_rate_filter_2.output.bed,
        bim=rules.snp_call_rate_filter_2.output.bim,
        fam=rules.snp_call_rate_filter_2.output.fam,
    params:
        maf="{maf}",
        out_prefix="sample_level/call_rate_2/samples_maf{maf}",
    output:
        bed=temp("sample_level/call_rate_2/samples_maf{maf}.bed"),
        bim=temp("sample_level/call_rate_2/samples_maf{maf}.bim"),
        fam=temp("sample_level/call_rate_2/samples_maf{maf}.fam"),
        nosex=temp("sample_level/call_rate_2/samples_maf{maf}.nosex"),
    log:
        "sample_level/call_rate_2/samples_maf{maf}.log",


use rule ld from plink as sample_level_ld_estimate with:
    input:
        bed=rules.sample_level_maf_filter.output.bed,
        bim=rules.sample_level_maf_filter.output.bim,
        fam=rules.sample_level_maf_filter.output.fam,
    params:
        r2="{ld}",  # r2 threshold: currently 0.1
        out_prefix="sample_level/call_rate_2/samples_maf{maf}_ld{ld}",
    output:
        to_keep=temp("sample_level/call_rate_2/samples_maf{maf}_ld{ld}.prune.in"),  # Markers in approx. linkage equilibrium
        to_remove=temp("sample_level/call_rate_2/samples_maf{maf}_ld{ld}.prune.out"),  # Markers in LD
        nosex=temp("sample_level/call_rate_2/samples_maf{maf}_ld{ld}.nosex"),
    log:
        "sample_level/call_rate_2/samples_maf{maf}_ld{ld}.log",


use rule ld_filter from plink as sample_level_ld_pruned with:
    input:
        bed=rules.sample_level_maf_filter.output.bed,
        bim=rules.sample_level_maf_filter.output.bim,
        fam=rules.sample_level_maf_filter.output.fam,
        to_keep=rules.sample_level_ld_estimate.output.to_keep,
    params:
        out_prefix="sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned",
    output:
        bed="sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.bed",
        bim="sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.bim",
        fam="sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.fam",
        nosex="sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.nosex",
    log:
        "sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.log",


use rule genome from plink as sample_level_ibd with:
    input:
        bed=rules.sample_level_ld_pruned.output.bed,
        bim=rules.sample_level_ld_pruned.output.bim,
        fam=rules.sample_level_ld_pruned.output.fam,
    params:
        ibd_min=cfg.config.software_params.ibd_pi_hat_min,
        ibd_max=cfg.config.software_params.ibd_pi_hat_max,
        out_prefix="sample_level/call_rate_2/samples_maf{maf}_ld{ld}",
    output:
        "sample_level/call_rate_2/samples_maf{maf}_ld{ld}.genome",


rule sample_concordance_plink:
    """Parse PLINK IBD file and calculate sample concordance.

    Where, sample concordance is `IBS2 / (IBS0 + IBS1 + IBS2)`.
    """
    input:
        rules.sample_level_ibd.output[0].format(
            maf=cfg.config.software_params.maf_for_ibd,
            ld=cfg.config.software_params.ld_prune_r2,
        ),
    params:
        concordance_threshold=cfg.config.software_params.dup_concordance_cutoff,
        pi_hat_threshold=cfg.config.software_params.pi_hat_threshold,
    output:
        "sample_level/concordance/plink.csv",
    script:
        "../scripts/concordance_table.py"


use rule extract_fingerprint_snps from graf as graf_extract_fingerprint_snps with:
    input:
        bed=rules.update_samples_to_1kg_rsIDs.output.bed,
        bim=rules.update_samples_to_1kg_rsIDs.output.bim,
        fam=rules.update_samples_to_1kg_rsIDs.output.fam,
    params:
        out_prefix="sample_level/call_rate_2/samples_1kg_rsID",
    output:
        "sample_level/call_rate_2/samples_1kg_rsID.fpg",
    log:
        "sample_level/call_rate_2/samples_1kg_rsID.fpg.log",


use rule relatedness from graf as sample_concordance_graf with:
    input:
        rules.graf_extract_fingerprint_snps.output[0],
    output:
        "sample_level/concordance/graf.tsv",
    log:
        "sample_level/concordance/graf.log",


rule sample_concordance_king:
    input:
        bed=rules.snp_call_rate_filter_2.output.bed,
        bim=rules.snp_call_rate_filter_2.output.bim,
        fam=rules.snp_call_rate_filter_2.output.fam,
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
        # king does not always output all files so touch for snakemake
        "king -b {input.bed} --related --degree 3 --prefix {params.out_prefix} --cpu {threads} > {log} 2>&1 "
        "&& touch {output.within_family} "
        "&& touch {output.between_family} "
        "&& touch {output.within_family_X} "
        "&& touch {output.between_family_X} "
        "&& touch {output.segments} "


rule sample_concordance_summary:
    input:
        sample_sheet_csv="cgr_sample_sheet.csv",
        plink_file=rules.sample_concordance_plink.output[0],
        graf_file=rules.sample_concordance_graf.output[0],
        king_file=rules.sample_concordance_king.output.between_family,
    output:
        "sample_level/concordance/summary.csv",
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * 2 * attempt,
        time_hr=lambda wildcards, attempt: 2 * attempt,
    script:
        "../scripts/sample_concordance.py"


rule split_sample_concordance:
    input:
        rules.sample_concordance_summary.output[0],
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
use rule populations from graf as graf_populations with:
    input:
        rules.graf_extract_fingerprint_snps.output[0],
    output:
        "sample_level/ancestry/graf_populations.txt",
    log:
        "sample_level/ancestry/graf_populations.log",


use rule ancestry from graf as graf_ancestry with:
    input:
        rules.graf_populations.output[0],
    output:
        "sample_level/ancestry/graf_ancestry.txt",


# -------------------------------------------------------------------------------
# SNP Summary Table
# -------------------------------------------------------------------------------
rule snp_qc_table:
    input:
        initial=rules.plink_call_rate_initial.output.lmiss,
        cr1=rules.plink_call_rate_post1.output.lmiss,
        cr2=rules.plink_call_rate_post2.output.lmiss,
        thousand_genomes=rules.update_samples_to_1kg_rsIDs.output.id_map,
    output:
        "sample_level/snp_qc.csv",
    script:
        "../scripts/snp_qc_table.py"


# -------------------------------------------------------------------------------
# Sample Summary Table
# -------------------------------------------------------------------------------
use rule sexcheck from plink as sample_level_sexcheck with:
    input:
        bed=rules.snp_call_rate_filter_1.output.bed,
        bim=rules.snp_call_rate_filter_1.output.bim,
        fam=rules.snp_call_rate_filter_1.output.fam,
    params:
        out_prefix="sample_level/call_rate_1/samples",
    output:
        "sample_level/call_rate_1/samples.sexcheck",


def _contam(wildcards):
    uf, wf = cfg.config.user_files, cfg.config.workflow_params
    if use_contamination:
        return rules.sample_contamination_verifyIDintensity.output[0]
    return []


def _intensity(wildcards):
    uf, wf = cfg.config.user_files, cfg.config.workflow_params
    if use_contamination:
        return "sample_level/contamination/median_idat_intensity.csv"
    return []


rule sample_qc_table:
    input:
        sample_sheet_csv="cgr_sample_sheet.csv",
        imiss_start=rules.plink_call_rate_initial.output.imiss,
        imiss_cr1=rules.plink_call_rate_post1.output.imiss,
        imiss_cr2=rules.plink_call_rate_post2.output.imiss,
        sexcheck_cr1=rules.sample_level_sexcheck.output[0],
        ancestry=rules.graf_ancestry.output[0],
        sample_concordance_csv=rules.sample_concordance_summary.output[0],
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
