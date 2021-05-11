from cgr_gwas_qc import load_config

cfg = load_config()
_targets = []


################################################################################
# Sub-workflows
################################################################################
subworkflow entry_points:
    snakefile:
        cfg.subworkflow("entry_points")
    workdir:
        cfg.root.as_posix()


_targets.append(entry_points("subworkflow_complete/entry_points"))


subworkflow contamination:
    snakefile:
        cfg.subworkflow("contamination")
    workdir:
        cfg.root.as_posix()


if (
    cfg.config.user_files.idat_pattern
    and cfg.config.user_files.gtc_pattern
    and cfg.config.workflow_params.remove_contam
):
    _targets.append(contamination("subworkflow_complete/contamination"))
    _targets.append("sample_level/contamination/summary.csv")


################################################################################
# Sample QC Targets
################################################################################
_targets.extend(
    [
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
)


rule all_sample_qc:
    input:
        _targets,
    output:
        touch("subworkflow_complete/sample_qc"),


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
    config:
        {}


module graf:
    snakefile:
        cfg.modules("graf")
    config:
        {}


use rule update_snps_to_1kg_rsID from thousand_genomes


################################################################################
# Workflow Rules
################################################################################
# -------------------------------------------------------------------------------
# Call Rate Filters
# -------------------------------------------------------------------------------
use rule sample_call_rate_filter from plink as sample_call_rate_filter_1 with:
    input:
        bed=entry_points("sample_level/samples.bed"),
        bim=entry_points("sample_level/samples.bim"),
        fam=entry_points("sample_level/samples.fam"),
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
    group:
        "call_rate_filters"


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
    group:
        "call_rate_filters"


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
    group:
        "call_rate_filters"


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
    group:
        "call_rate_filters"


use rule miss from plink as plink_miss_initial with:
    input:
        bed=entry_points("sample_level/samples.bed"),
        bim=entry_points("sample_level/samples.bim"),
        fam=entry_points("sample_level/samples.fam"),
    params:
        out_prefix="sample_level/samples",
    output:
        imiss="sample_level/samples.imiss",
        lmiss="sample_level/samples.lmiss",
    group:
        "call_rate_filters"


use rule miss from plink as plink_miss_cr with:
    input:
        bed="sample_level/call_rate_{cr}/samples.bed",
        bim="sample_level/call_rate_{cr}/samples.bim",
        fam="sample_level/call_rate_{cr}/samples.fam",
    params:
        out_prefix="sample_level/call_rate_{cr}/samples",
    output:
        imiss="sample_level/call_rate_{cr}/samples.imiss",
        lmiss="sample_level/call_rate_{cr}/samples.lmiss",
    wildcard_constraints:
        cr="1|2",
    group:
        "call_rate_filters"


# -------------------------------------------------------------------------------
# Contamination
# -------------------------------------------------------------------------------
rule sample_contamination_verifyIDintensity:
    """Aggregate sample contamination scores.

    Aggregates sample level contamination scores into a single file (each
    row is a sample). The script sets ``%Mix`` to ``NA`` if the intensity
    is below the threshold and the file is not in the ``imiss`` file.
    """
    input:
        contamination_file=contamination("sample_level/contamination/verifyIDintensity.csv"),
        median_intensity_file=contamination(
            "sample_level/contamination/median_idat_intensity.csv"
        ),
        imiss_file="sample_level/call_rate_2/samples.imiss",
    params:
        intensity_threshold=cfg.config.software_params.intensity_threshold,
        contam_threshold=cfg.config.software_params.contam_threshold,
    output:
        "sample_level/contamination/summary.csv",
    group:
        "sample_qc"
    script:
        "../scripts/agg_contamination.py"


# -------------------------------------------------------------------------------
# Replicate Concordance
# -------------------------------------------------------------------------------
use rule maf_filter from plink as plink_* with:
    group:
        "replicate_concordance"


use rule ld from plink as plink_* with:
    group:
        "replicate_concordance"


use rule ld_filter from plink as plink_* with:
    group:
        "replicate_concordance"


use rule keep_bfile from plink as plink_* with:
    group:
        "replicate_concordance"


use rule genome from plink as plink_* with:
    group:
        "replicate_concordance"


rule sample_concordance_plink:
    """Parse PLINK IBD file and calculate sample concordance.

    Where, sample concordance is `IBS2 / (IBS0 + IBS1 + IBS2)`.
    """
    input:
        "sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.keep.genome".format(
            maf=cfg.config.software_params.maf_for_ibd,
            ld=cfg.config.software_params.ld_prune_r2,
        ),
    params:
        concordance_threshold=cfg.config.software_params.dup_concordance_cutoff,
        pi_hat_threshold=cfg.config.software_params.pi_hat_threshold,
    output:
        "sample_level/concordance/plink.csv",
    group:
        "replicate_concordance"
    script:
        "../scripts/concordance_table.py"


use rule extract_fingerprint_snps from graf as graf_extract_fpg with:
    group:
        "replicate_concordance"


use rule relatedness from graf as sample_concordance_graf with:
    input:
        "sample_level/call_rate_2/samples_1kg_rsID.fpg",
    output:
        "sample_level/concordance/graf.tsv",
    log:
        "sample_level/concordance/graf.log",
    group:
        "replicate_concordance"


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
    group:
        "replicate_concordance"
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
    group:
        "replicate_concordance"
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
    group:
        "replicate_concordance"
    script:
        "../scripts/split_sample_concordance.py"


# -------------------------------------------------------------------------------
# Ancestry
# -------------------------------------------------------------------------------
use rule populations from graf as graf_populations with:
    input:
        "sample_level/call_rate_2/samples_1kg_rsID.fpg",
    output:
        "sample_level/ancestry/graf_populations.txt",
    log:
        "sample_level/ancestry/graf_populations.log",
    group:
        "ancestry"


use rule ancestry from graf as graf_ancestry with:
    input:
        rules.graf_populations.output[0],
    output:
        "sample_level/ancestry/graf_ancestry.txt",
    group:
        "ancestry"


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
    group:
        "sample_qc"
    script:
        "../scripts/snp_qc_table.py"


# -------------------------------------------------------------------------------
# Sample Summary Table
# -------------------------------------------------------------------------------
use rule sexcheck from plink as plink_sexcheck with:
    group:
        "sample_qc"


def _contam(wildcards):
    uf, wf = cfg.config.user_files, cfg.config.workflow_params
    if uf.idat_pattern and uf.gtc_pattern and wf.remove_contam:
        return rules.sample_contamination_verifyIDintensity.output[0]
    return []


def _intensity(wildcards):
    uf, wf = cfg.config.user_files, cfg.config.workflow_params
    if uf.idat_pattern and uf.gtc_pattern and wf.remove_contam:
        return contamination("sample_level/contamination/median_idat_intensity.csv")
    return []


rule sample_qc_table:
    input:
        sample_sheet_csv="cgr_sample_sheet.csv",
        imiss_start="sample_level/samples.imiss",
        imiss_cr1="sample_level/call_rate_1/samples.imiss",
        imiss_cr2="sample_level/call_rate_2/samples.imiss",
        sexcheck_cr1="sample_level/call_rate_1/samples.sexcheck",
        ancestry=rules.graf_ancestry.output[0],
        sample_concordance_csv=rules.sample_concordance_summary.output[0],
        contam=_contam,
        intensity=_intensity,
    params:
        remove_contam=cfg.config.workflow_params.remove_contam,
        remove_rep_discordant=cfg.config.workflow_params.remove_rep_discordant,
    output:
        "sample_level/sample_qc.csv",
    group:
        "sample_qc"
    script:
        "../scripts/sample_qc_table.py"


rule sample_qc_summary_stats:
    input:
        rules.sample_qc_table.output[0],
    output:
        "sample_level/summary_stats.txt",
    group:
        "sample_qc"
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
    group:
        "sample_qc"
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
    group:
        "sample_qc"
    script:
        "../scripts/plot_call_rate.py"


rule plot_chrx_inbreeding:
    input:
        rules.sample_qc_table.output[0],
    output:
        "sample_level/chrx_inbreeding.png",
    group:
        "sample_qc"
    script:
        "../scripts/plot_chrx_inbreeding.py"


rule plot_ancestry:
    input:
        rules.sample_qc_table.output[0],
    output:
        "sample_level/ancestry.png",
    group:
        "sample_qc"
    script:
        "../scripts/plot_ancestry.py"
