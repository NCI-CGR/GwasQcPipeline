from more_itertools import flatten

from cgr_gwas_qc import load_config
from cgr_gwas_qc.workflow.scripts import subject_qc_table


cfg = load_config()


wildcard_constraints:
    population="[\w_]+",


################################################################################
# Subject QC Targets
################################################################################
targets = [
    "subject_level/concordance.csv",
    "subject_level/population_qc.csv",
    "subject_level/.population_plots.done",
    "subject_level/.control_plots.done",
]


rule all_subject_qc:
    input:
        targets,


################################################################################
# Imports
################################################################################
module plink:
    snakefile:
        cfg.modules("plink")
    config:
        {}


module eigensoft:
    snakefile:
        cfg.modules("eigensoft")
    config:
        {}


################################################################################
# Workflow Rules
################################################################################
# -------------------------------------------------------------------------------
# Subject Level Analysis
# -------------------------------------------------------------------------------
# Subject Summary Table
rule subject_qc_table:
    input:
        sample_qc_csv="sample_level/sample_qc.csv",
        sample_concordance_csv="sample_level/concordance/summary.csv",
    params:
        remove_sex_discordant=cfg.config.workflow_params.remove_sex_discordant,
        remove_unexpected_rep=cfg.config.workflow_params.remove_unexpected_rep,
    output:
        "subject_level/subject_qc.csv",
    group:
        "select_subjects"
    script:
        "../scripts/subject_qc_table.py"


checkpoint select_subjects_for_analysis:
    """Pull out subjects that have no analytic exclusions."""
    input:
        rules.subject_qc_table.output[0],
    output:
        "subject_level/subjects_for_analysis.csv",
    group:
        "select_subjects"
    run:
        (
            subject_qc_table.read(input[0])
            .query("not subject_analytic_exclusion")
            .to_csv(output[0], index=False)
        )


# Convert Sample Level data to Subject Level
rule selected_Subject_IDs:
    """Create mapping file for plink to convert Sample_IDs to Subject_IDs."""
    input:
        rules.select_subjects_for_analysis.output[0],
    output:
        selected=temp("subject_level/selected_Subject_IDs.txt"),
        sample2subject=temp("subject_level/sample_to_subject.txt"),
    group:
        "select_subjects"
    run:
        qc = subject_qc_table.read(input[0])
        qc.reindex(["Sample_ID", "Sample_ID"], axis=1).to_csv(
            output.selected, sep=" ", index=False, header=False
        )
        qc.reindex(
            ["Sample_ID", "Sample_ID", "Group_By_Subject_ID", "Group_By_Subject_ID"],
            axis=1,
        ).to_csv(output.sample2subject, sep=" ", index=False, header=False)


use rule keep_ids from plink as pull_selected_subjects with:
    input:
        bed="sample_level/call_rate_2/samples.bed",
        bim="sample_level/call_rate_2/samples.bim",
        fam="sample_level/call_rate_2/samples.fam",
        to_keep=rules.selected_Subject_IDs.output.selected,
    params:
        out_prefix="subject_level/samples",
    output:
        bed="subject_level/samples.bed",
        bim="subject_level/samples.bim",
        fam="subject_level/samples.fam",
        nosex="subject_level/samples.nosex",
    log:
        "subject_level/samples.log",
    group:
        "select_subjects"


use rule rename_ids from plink as rename_Sample_ID_to_Subject_ID with:
    input:
        bed=rules.pull_selected_subjects.output.bed,
        bim=rules.pull_selected_subjects.output.bim,
        fam=rules.pull_selected_subjects.output.fam,
        id_map=rules.selected_Subject_IDs.output.sample2subject,
    params:
        out_prefix="subject_level/subjects",
    output:
        bed="subject_level/subjects.bed",
        bim="subject_level/subjects.bim",
        fam="subject_level/subjects.fam",
        nosex="subject_level/subjects.nosex",
    log:
        "subject_level/subjects.log",
    group:
        "select_subjects"


# -------------------------------------------------------------------------------
# Population Level Analysis
# -------------------------------------------------------------------------------
def _get_populations(wildcards):
    """Use checkpoints to pull out a list of populations"""
    filename = checkpoints.select_subjects_for_analysis.get(**wildcards).output[0]
    min_num_subjects = cfg.config.workflow_params.minimum_pop_subjects
    return (
        subject_qc_table.read(filename)
        .groupby("Ancestry")
        .size()
        .pipe(lambda x: x[x >= min_num_subjects])  # defaults to 50
        .index.tolist()
    )


# Split Subjects by Populations
rule population_Subject_IDs:
    input:
        rules.select_subjects_for_analysis.output[0],
    output:
        "subject_level/{population}/subjects.txt",
    group:
        "{population}"
    run:
        (
            subject_qc_table.read(input[0])
            .query("Ancestry == @wildcards.population")
            .reindex(["Group_By_Subject_ID", "Group_By_Subject_ID"], axis=1)
            .to_csv(output[0], sep=" ", index=False, header=False)
        )


# Within Population Subject Concordance
use rule keep_ids from plink as pull_population_level_subjects with:
    input:
        bed=rules.rename_Sample_ID_to_Subject_ID.output.bed,
        bim=rules.rename_Sample_ID_to_Subject_ID.output.bim,
        fam=rules.rename_Sample_ID_to_Subject_ID.output.fam,
        to_keep=rules.population_Subject_IDs.output[0],
    params:
        out_prefix="subject_level/{population}/subjects",
    output:
        bed="subject_level/{population}/subjects.bed",
        bim="subject_level/{population}/subjects.bim",
        fam="subject_level/{population}/subjects.fam",
        nosex="subject_level/{population}/subjects.nosex",
    log:
        "subject_level/{population}/subjects.log",
    group:
        "{population}"


use rule maf_filter from plink as population_level_maf_filter with:
    input:
        bed=rules.pull_population_level_subjects.output.bed,
        bim=rules.pull_population_level_subjects.output.bim,
        fam=rules.pull_population_level_subjects.output.fam,
    params:
        maf="{maf}",
        out_prefix="subject_level/{population}/subjects_maf{maf}",
    output:
        bed=temp("subject_level/{population}/subjects_maf{maf}.bed"),
        bim=temp("subject_level/{population}/subjects_maf{maf}.bim"),
        fam=temp("subject_level/{population}/subjects_maf{maf}.fam"),
        nosex=temp("subject_level/{population}/subjects_maf{maf}.nosex"),
    log:
        "subject_level/{population}/subjects_maf{maf}.log",
    group:
        "{population}"


use rule ld from plink as population_level_ld_estimate with:
    input:
        bed=rules.population_level_maf_filter.output.bed,
        bim=rules.population_level_maf_filter.output.bim,
        fam=rules.population_level_maf_filter.output.fam,
    params:
        r2="{ld}",  # r2 threshold: currently 0.1
        out_prefix="subject_level/{population}/subjects_maf{maf}_ld{ld}",
    output:
        # Markers in approx. linkage equilibrium
        to_keep=temp("subject_level/{population}/subjects_maf{maf}_ld{ld}.prune.in"),
        # Markers in LD
        to_remove=temp("subject_level/{population}/subjects_maf{maf}_ld{ld}.prune.out"),
        nosex=temp("subject_level/{population}/subjects_maf{maf}_ld{ld}.nosex"),
    log:
        "subject_level/{population}/subjects_maf{maf}_ld{ld}.log",
    group:
        "{population}"


use rule ld_filter from plink as population_level_ld_pruned with:
    input:
        bed=rules.population_level_maf_filter.output.bed,
        bim=rules.population_level_maf_filter.output.bim,
        fam=rules.population_level_maf_filter.output.fam,
        to_keep=rules.population_level_ld_estimate.output.to_keep,
    params:
        out_prefix="subject_level/{population}/subjects_maf{maf}_ld{ld}_pruned",
    output:
        bed="subject_level/{population}/subjects_maf{maf}_ld{ld}_pruned.bed",
        bim="subject_level/{population}/subjects_maf{maf}_ld{ld}_pruned.bim",
        fam="subject_level/{population}/subjects_maf{maf}_ld{ld}_pruned.fam",
        nosex="subject_level/{population}/subjects_maf{maf}_ld{ld}_pruned.nosex",
    log:
        "subject_level/{population}/subjects_maf{maf}_ld{ld}_pruned.log",
    group:
        "{population}"


use rule genome from plink as population_level_ibd with:
    input:
        bed=rules.population_level_ld_pruned.output.bed,
        bim=rules.population_level_ld_pruned.output.bim,
        fam=rules.population_level_ld_pruned.output.fam,
    params:
        ibd_min=cfg.config.software_params.ibd_pi_hat_min,
        ibd_max=cfg.config.software_params.ibd_pi_hat_max,
        out_prefix="subject_level/{population}/subjects_maf{maf}_ld{ld}",
    output:
        "subject_level/{population}/subjects_maf{maf}_ld{ld}.genome",
    group:
        "{population}"


rule population_level_concordance_plink:
    input:
        rules.population_level_ibd.output[0],
    params:
        concordance_threshold=cfg.config.software_params.dup_concordance_cutoff,
        pi_hat_threshold=cfg.config.software_params.pi_hat_threshold,
    output:
        "subject_level/{population}/subjects_maf{maf}_ld{ld}.concordance.csv",
    group:
        "{population}"
    script:
        "../scripts/concordance_table.py"


def _population_concordance_files(wildcards):
    populations = _get_populations(wildcards)

    if not populations:
        return []

    return expand(
        rules.population_level_concordance_plink.output[0],
        population=populations,
        maf=cfg.config.software_params.maf_for_ibd,
        ld=cfg.config.software_params.ld_prune_r2,
    )


rule agg_population_concordance:
    input:
        sample_sheet_csv="cgr_sample_sheet.csv",
        ibd_files=_population_concordance_files,
    output:
        "subject_level/concordance.csv",
    group:
        "subject_qc"
    script:
        "../scripts/agg_population_concordance.py"


# Remove Related Subjects
rule population_level_related_subjects:
    input:
        expand(
            rules.population_level_concordance_plink.output[0],
            maf=cfg.config.software_params.maf_for_ibd,
            ld=cfg.config.software_params.ld_prune_r2,
            allow_missing=True,
        ),
    output:
        relatives="subject_level/{population}/relatives.csv",
        to_remove="subject_level/{population}/related_subjects_to_remove.txt",
    group:
        "{population}"
    script:
        "../scripts/related_subjects.py"


use rule remove_ids from plink as population_level_remove_related_subjects with:
    input:
        bed=rules.pull_population_level_subjects.output.bed,
        bim=rules.pull_population_level_subjects.output.bim,
        fam=rules.pull_population_level_subjects.output.fam,
        to_remove=rules.population_level_related_subjects.output.to_remove,
    params:
        out_prefix="subject_level/{population}/subjects_unrelated",
    output:
        bed="subject_level/{population}/subjects_unrelated.bed",
        bim="subject_level/{population}/subjects_unrelated.bim",
        fam="subject_level/{population}/subjects_unrelated.fam",
        nosex="subject_level/{population}/subjects_unrelated.nosex",
    log:
        "subject_level/{population}/subjects_unrelated.log",
    group:
        "{population}"


# PCA
use rule maf_filter from plink as population_level_unrelated_maf_filter with:
    input:
        bed=rules.population_level_remove_related_subjects.output.bed,
        bim=rules.population_level_remove_related_subjects.output.bim,
        fam=rules.population_level_remove_related_subjects.output.fam,
    params:
        maf="{maf}",
        out_prefix="subject_level/{population}/subjects_unrelated_maf{maf}",
    output:
        bed=temp("subject_level/{population}/subjects_unrelated_maf{maf}.bed"),
        bim=temp("subject_level/{population}/subjects_unrelated_maf{maf}.bim"),
        fam=temp("subject_level/{population}/subjects_unrelated_maf{maf}.fam"),
        nosex=temp("subject_level/{population}/subjects_unrelated_maf{maf}.nosex"),
    log:
        "subject_level/{population}/subjects_unrelated_maf{maf}.log",
    group:
        "{population}"


use rule ld from plink as population_level_unrelated_ld_estimate with:
    input:
        bed=rules.population_level_unrelated_maf_filter.output.bed,
        bim=rules.population_level_unrelated_maf_filter.output.bim,
        fam=rules.population_level_unrelated_maf_filter.output.fam,
    params:
        r2="{ld}",  # r2 threshold: currently 0.1
        out_prefix="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}",
    output:
        # Markers in approx. linkage equilibrium
        to_keep=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.prune.in"),
        # Markers in LD
        to_remove=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.prune.out"),
        nosex="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.nosex",
    log:
        "subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.log",
    group:
        "{population}"


use rule ld_filter from plink as population_level_unrelated_ld_pruned with:
    input:
        bed=rules.population_level_unrelated_maf_filter.output.bed,
        bim=rules.population_level_unrelated_maf_filter.output.bim,
        fam=rules.population_level_unrelated_maf_filter.output.fam,
        to_keep=rules.population_level_unrelated_ld_estimate.output.to_keep,
    params:
        out_prefix="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_pruned",
    output:
        bed="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_pruned.bed",
        bim="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_pruned.bim",
        fam="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_pruned.fam",
        nosex="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_pruned.nosex",
    log:
        "subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_pruned.log",
    group:
        "{population}"


use rule bed_to_ped from plink as population_level_unrelated_bed_to_ped with:
    input:
        bed=rules.population_level_unrelated_ld_pruned.output.bed,
        bim=rules.population_level_unrelated_ld_pruned.output.bim,
        fam=rules.population_level_unrelated_ld_pruned.output.fam,
    output:
        ped=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.ped"),
        map_=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.map"),
    params:
        out_prefix="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}",
    log:
        "subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.log",
    group:
        "{population}"


use rule convert from eigensoft as population_level_unrelated_convert_to_eigensoft with:
    input:
        ped=rules.population_level_unrelated_bed_to_ped.output.ped,
        map_=rules.population_level_unrelated_bed_to_ped.output.map_,
    output:
        par=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.convert.par"),
        gen=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.gen"),
        snp=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.snp"),
        ind=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.ind"),
    group:
        "{population}"


use rule smartpca from eigensoft as population_level_unrelated_smartpca with:
    input:
        gen=rules.population_level_unrelated_convert_to_eigensoft.output.gen,
        snp=rules.population_level_unrelated_convert_to_eigensoft.output.snp,
        ind=rules.population_level_unrelated_convert_to_eigensoft.output.ind,
    output:
        par=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.pca.par"),
        eigenvec="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.eigenvec",
    group:
        "{population}"


rule plot_pca:
    input:
        qc_table=rules.select_subjects_for_analysis.output[0],
        eigenvec=expand(
            rules.population_level_unrelated_smartpca.output.eigenvec,
            maf=cfg.config.software_params.maf_for_ibd,
            ld=cfg.config.software_params.ld_prune_r2,
            allow_missing=True,
        )[0],
    params:
        population="{population}",
    output:
        "subject_level/pca_plots/{population}.png",
    group:
        "{population}"
    script:
        "../scripts/plot_pca.py"


# Autosomal Heterozygosity (see plink_stats.smk 'plink_stats_het')
use rule het from plink as population_level_autosomal_heterozygosity with:
    input:
        bed=rules.pull_population_level_subjects.output.bed,
        bim=rules.pull_population_level_subjects.output.bim,
        fam=rules.pull_population_level_subjects.output.fam,
    params:
        out_prefix="subject_level/{population}/subjects",
    output:
        "subject_level/{population}/subjects.het",
    group:
        "{population}"


rule plot_autosomal_heterozygosity:
    input:
        qc_table=rules.select_subjects_for_analysis.output[0],
        het=rules.population_level_autosomal_heterozygosity.output[0],
    params:
        population="{population}",
        threshold=cfg.config.software_params.autosomal_het_threshold,
    output:
        "subject_level/autosomal_heterozygosity_plots/{population}.png",
    group:
        "{population}"
    script:
        "../scripts/plot_autosomal_heterozygosity.py"


# Population Summary Table
rule per_population_qc_table:
    input:
        relatives=rules.population_level_related_subjects.output.relatives,
        pca=expand(
            rules.population_level_unrelated_smartpca.output.eigenvec,
            maf=cfg.config.software_params.maf_for_ibd,
            ld=cfg.config.software_params.ld_prune_r2,
            allow_missing=True,
        )[0],
        autosomal_het=rules.population_level_autosomal_heterozygosity.output[0],
    params:
        population="{population}",
        threshold=cfg.config.software_params.autosomal_het_threshold,
    output:
        "subject_level/{population}/qc.csv",
    group:
        "{population}"
    script:
        "../scripts/population_qc_table.py"


def _population_qc_tables(wildcards):
    populations = _get_populations(wildcards)

    if not populations:
        return []

    return expand(rules.per_population_qc_table.output[0], population=populations)


rule agg_population_qc_tables:
    input:
        subject_qc_table=rules.select_subjects_for_analysis.output[0],
        population_qc_tables=_population_qc_tables,
    output:
        "subject_level/population_qc.csv",
    group:
        "subject_qc"
    script:
        "../scripts/agg_population_qc_tables.py"


def _population_plots(wildcards):
    populations = _get_populations(wildcards)

    if not populations:
        return []

    return flatten(
        [
            expand(rules.plot_pca.output[0], population=populations),
            expand(rules.plot_autosomal_heterozygosity.output[0], population=populations),
        ]
    )


rule agg_population_plots:
    input:
        _population_plots,
    output:
        touch("subject_level/.population_plots.done"),


# -------------------------------------------------------------------------------
# Control Analysis (HWE)
# -------------------------------------------------------------------------------
rule population_controls_Subject_IDs:
    input:
        rules.select_subjects_for_analysis.output[0],
    output:
        "subject_level/{population}/controls.txt",
    group:
        "{population}"
    run:
        (
            subject_qc_table.read(input[0])
            .query("Ancestry == @wildcards.population & case_control == 'Control'")
            .reindex(["Group_By_Subject_ID", "Group_By_Subject_ID"], axis=1)
            .to_csv(output[0], sep=" ", index=False, header=False)
        )


use rule keep_ids from plink as pull_population_unrelated_controls with:
    input:
        bed=rules.population_level_remove_related_subjects.output.bed,
        bim=rules.population_level_remove_related_subjects.output.bim,
        fam=rules.population_level_remove_related_subjects.output.fam,
        to_keep=rules.population_controls_Subject_IDs.output[0],
    params:
        out_prefix="subject_level/{population}/controls_unrelated",
    output:
        bed="subject_level/{population}/controls_unrelated.bed",
        bim="subject_level/{population}/controls_unrelated.bim",
        fam="subject_level/{population}/controls_unrelated.fam",
        nosex="subject_level/{population}/controls_unrelated.nosex",
    group:
        "{population}"
    log:
        "subject_level/{population}/controls_unrelated.log",


use rule maf_filter from plink as population_level_unrelated_controls_maf_filter with:
    input:
        bed=rules.pull_population_unrelated_controls.output.bed,
        bim=rules.pull_population_unrelated_controls.output.bim,
        fam=rules.pull_population_unrelated_controls.output.fam,
    params:
        maf="{maf}",
        out_prefix="subject_level/{population}/controls_unrelated_maf{maf}",
    output:
        bed=temp("subject_level/{population}/controls_unrelated_maf{maf}.bed"),
        bim=temp("subject_level/{population}/controls_unrelated_maf{maf}.bim"),
        fam=temp("subject_level/{population}/controls_unrelated_maf{maf}.fam"),
        nosex=temp("subject_level/{population}/controls_unrelated_maf{maf}.nosex"),
    log:
        "subject_level/{population}/controls_unrelated_maf{maf}.log",
    group:
        "{population}"


use rule snps_only_filter from plink as keep_only_snps with:
    input:
        bed=rules.population_level_unrelated_controls_maf_filter.output.bed,
        bim=rules.population_level_unrelated_controls_maf_filter.output.bim,
        fam=rules.population_level_unrelated_controls_maf_filter.output.fam,
    params:
        out_prefix="subject_level/{population}/controls_unrelated_maf{maf}_snps",
    output:
        bed=temp("subject_level/{population}/controls_unrelated_maf{maf}_snps.bed"),
        bim=temp("subject_level/{population}/controls_unrelated_maf{maf}_snps.bim"),
        fam=temp("subject_level/{population}/controls_unrelated_maf{maf}_snps.fam"),
        nosex=temp("subject_level/{population}/controls_unrelated_maf{maf}_snps.nosex"),
    log:
        "subject_level/{population}/controls_unrelated_maf{maf}_snps.log",
    group:
        "{population}"


use rule autosome_only_filter from plink as keep_only_autosomal_variants with:
    input:
        bed=rules.keep_only_snps.output.bed,
        bim=rules.keep_only_snps.output.bim,
        fam=rules.keep_only_snps.output.fam,
    params:
        out_prefix="subject_level/{population}/controls_unrelated_maf{maf}_snps_autosome",
    output:
        bed=temp("subject_level/{population}/controls_unrelated_maf{maf}_snps_autosome.bed"),
        bim=temp("subject_level/{population}/controls_unrelated_maf{maf}_snps_autosome.bim"),
        fam=temp("subject_level/{population}/controls_unrelated_maf{maf}_snps_autosome.fam"),
        nosex=temp("subject_level/{population}/controls_unrelated_maf{maf}_snps_autosome.nosex"),
    log:
        "subject_level/{population}/controls_unrelated_maf{maf}_snps_autosome.log",
    group:
        "{population}"


use rule hwe from plink as population_level_unrelated_controls_hwe with:
    input:
        bed=rules.keep_only_autosomal_variants.output.bed,
        bim=rules.keep_only_autosomal_variants.output.bim,
        fam=rules.keep_only_autosomal_variants.output.fam,
    params:
        out_prefix="subject_level/{population}/controls_unrelated_maf{maf}_snps_autosome",
    output:
        "subject_level/{population}/controls_unrelated_maf{maf}_snps_autosome.hwe",
    group:
        "{population}"


rule plot_hwe:
    input:
        expand(
            rules.population_level_unrelated_controls_hwe.output[0],
            maf=cfg.config.software_params.maf_for_hwe,
            allow_missing=True,
        ),
    params:
        population="{population}",
    output:
        "subject_level/hwe_plots/{population}.png",
    group:
        "{population}"
    script:
        "../scripts/plot_hwe.py"


def _control_plots(wildcards):
    filename = checkpoints.select_subjects_for_analysis.get(**wildcards).output[0]
    min_num_subjects = cfg.config.workflow_params.control_hwp_threshold
    populations = (
        subject_qc_table.read(filename)
        .query("case_control == 'Control'")
        .groupby("Ancestry")
        .size()
        .pipe(lambda x: x[x >= min_num_subjects])  # defaults to 50
        .index.tolist()
    )

    if not populations:
        return []

    return expand(rules.plot_hwe.output[0], population=populations)


rule agg_control_plots:
    input:
        _control_plots,
    output:
        touch("subject_level/.control_plots.done"),
    group:
        "subject_qc"
