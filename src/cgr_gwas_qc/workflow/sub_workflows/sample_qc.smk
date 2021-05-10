from cgr_gwas_qc import load_config

cfg = load_config()


################################################################################
# All Targets
################################################################################
def _contamination_outputs(wildcards):
    """Only build contamination targets if we have IDAT & GTC & remove contam set"""
    if (
        cfg.config.user_files.idat_pattern
        and cfg.config.user_files.gtc_pattern
        and cfg.config.workflow_params.remove_contam
    ):
        return [
            "sample_level/median_idat_intensity.csv",
            "sample_level/contamination.csv",
        ]

    return []


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
        _contamination_outputs,


################################################################################
# Imports
################################################################################
include: cfg.modules("thousand_genomes")
include: cfg.subworkflow("entry_points")


module plink:
    snakefile:
        cfg.modules("plink")


use rule maf_filter, ld_filter, keep_bfile, miss, genome, ld, sexcheck from plink as plink_*


module graf:
    snakefile:
        cfg.modules("graf")


use rule extract_fingerprint_snps from graf


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
# Contamination
# -------------------------------------------------------------------------------
rule per_sample_median_idat_intensity:
    """Calculate median intensity of Red and Green channels."""
    input:
        red=lambda wc: cfg.expand(
            cfg.config.user_files.idat_pattern.red,
            query=f"Sample_ID == '{wc.Sample_ID}'",
        ),
        green=lambda wc: cfg.expand(
            cfg.config.user_files.idat_pattern.green,
            query=f"Sample_ID == '{wc.Sample_ID}'",
        ),
    output:
        temp(
            "sample_level/per_sample_median_idat_intensity/{Sample_ID}.{SentrixBarcode_A}.{SentrixPosition_A}.csv"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    group:
        "per_sample_median_idat_intensity"
    conda:
        cfg.conda("illuminaio")
    script:
        "../scripts/median_idat_intensity.R"


rule agg_median_idat_intensity:
    input:
        cfg.expand(rules.per_sample_median_idat_intensity.output[0]),
    output:
        "sample_level/median_idat_intensity.csv",
    run:
        import pandas as pd

        pd.concat([pd.read_csv(file_name) for file_name in input]).to_csv(
            output[0], index=False
        )


rule per_sample_gtc_to_adpc:
    """Converts a sample's GTC/BPM to an Illumina ADPC.BIN.

    This is the format required by ``verifyIDintensity``. The script also
    runs some sanity checks (intensities and normalized intensities > 0;
    genotypes are one of {0, 1, 2, 3}) while processing each file.

    .. warning::
        This is a submission hot spot creating 1 job per sample.
    """
    input:
        gtc_file=lambda wc: cfg.expand(
            cfg.config.user_files.gtc_pattern,
            query=f"Sample_ID == '{wc.Sample_ID}'",
        )[0],
        bpm_file=cfg.config.reference_files.illumina_manifest_file,
    output:
        temp("sample_level/per_sample_adpc/{Sample_ID}.adpc.bin"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    group:
        "per_sample_gtc_to_adpc"
    script:
        "../scripts/gtc2adpc.py"


rule per_sample_verifyIDintensity_contamination:
    """Find contaminated samples using allele intensities.

    Uses ``verifyIDintensity`` to find samples with allele intensities that deviate from the
    population.

    .. warning::
        This is a submission hot spot creating 1 job per sample.

    .. note::
        Here we are running ``verifyIDintensity`` in single sample mode. This software also has a
        multi-sample mode which may be faster and give better estimates. The problem with
        multi-sample mode is that it only works when you have a "large" number of samples.
    """
    input:
        adpc=rules.per_sample_gtc_to_adpc.output[0],
        abf=rules.pull_1KG_allele_b_freq.output.abf_file,
    params:
        snps=cfg.config.num_snps,
    output:
        temp("sample_level/per_sample_contamination_test/{Sample_ID}.contam.out"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    group:
        "per_sample_verifyIDintensity_contamination"
    conda:
        cfg.conda("verifyidintensity")
    shell:
        "verifyIDintensity -m {params.snps} -n 1 -b {input.abf} -v -p -i {input.adpc} > {output}"


rule agg_verifyIDintensity_contamination:
    """Aggregate sample contamination scores.

    Aggregates sample level contamination scores into a single file (each
    row is a sample). The script sets ``%Mix`` to ``NA`` if the intensity
    is below the threshold and the file is not in the ``imiss`` file.
    """
    input:
        contamination_files=cfg.expand(
            rules.per_sample_verifyIDintensity_contamination.output
        ),
        median_intensity_file=rules.agg_median_idat_intensity.output[0],
        imiss_file="sample_level/call_rate_2/samples.imiss",
    params:
        intensity_threshold=cfg.config.software_params.intensity_threshold,
        contam_threshold=cfg.config.software_params.contam_threshold,
    output:
        "sample_level/contamination.csv",
    script:
        "../scripts/agg_contamination.py"


# -------------------------------------------------------------------------------
# Replicate Concordance
# -------------------------------------------------------------------------------
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
    script:
        "../scripts/concordance_table.py"


use rule relatedness from graf as sample_concordance_graf with:
    input:
        "sample_level/call_rate_2/samples_1kg_rsID.fpg",
    output:
        "sample_level/concordance/graf.tsv",
    log:
        "sample_level/concordance/graf.log",


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
        "&& touch {output.segments} "
        # king does not always output all files so touch for snakemake


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
        "sample_level/call_rate_2/samples_1kg_rsID.fpg",
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
        return rules.agg_verifyIDintensity_contamination.output[0]
    return []


def _intensity(wildcards):
    uf, wf = cfg.config.user_files, cfg.config.workflow_params
    if uf.idat_pattern and uf.gtc_pattern and wf.remove_contam:
        return rules.agg_median_idat_intensity.output[0]
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
