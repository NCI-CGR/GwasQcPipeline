from more_itertools import flatten

from cgr_gwas_qc import load_config
from cgr_gwas_qc.workflow.scripts import subject_qc_table


cfg = load_config()


localrules:
    all_subject_qc,


wildcard_constraints:
    population="[\w_]+",


################################################################################
# All Targets
################################################################################
rule all_subject_qc:
    input:
        "subject_level/concordance.csv",
        "subject_level/population_qc.csv",
        "subject_level/.population_plots.done",
        "subject_level/.control_plots.done",


################################################################################
# Imports
################################################################################
subworkflow sample_qc:
    snakefile:
        cfg.subworkflow("sample_qc")
    workdir:
        cfg.root.as_posix()


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
        sample_qc_csv=sample_qc("sample_level/sample_qc.csv"),
        sample_concordance_csv=sample_qc("sample_level/concordance/summary.csv"),
    params:
        remove_sex_discordant=cfg.config.workflow_params.remove_sex_discordant,
        remove_unexpected_rep=cfg.config.workflow_params.remove_unexpected_rep,
    output:
        "subject_level/subject_qc.csv",
    group:
        "select_subjects"
    script:
        "../scripts/subject_qc_table.py"


checkpoint selected_subjects:
    """Pull out subjects that have no analytic exclusions."""
    input:
        rules.subject_qc_table.output[0],
    output:
        "subject_level/selected_subjects.csv",
    group:
        "select_subjects"
    run:
        (
            subject_qc_table.read(input[0])
            .query("not subject_analytic_exclusion")
            .to_csv(output[0], index=False)
        )


# Convert Sample Level data to Subject Level
rule list_of_subjects_to_analyze:
    """Create mapping file for plink to convert Sample_IDs to Subject_IDs."""
    input:
        rules.selected_subjects.output[0],
    output:
        selected=temp("subject_level/selected_subjects.txt"),
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


use rule keep_ids from plink as plink_subject_filter with:
    input:
        bed=sample_qc("sample_level/call_rate_2/samples.bed"),
        bim=sample_qc("sample_level/call_rate_2/samples.bim"),
        fam=sample_qc("sample_level/call_rate_2/samples.fam"),
        to_keep=rules.list_of_subjects_to_analyze.output.selected,
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


use rule rename_ids from plink as plink_convert_sample_id_to_subject_id with:
    input:
        bed=rules.plink_subject_filter.output.bed,
        bim=rules.plink_subject_filter.output.bim,
        fam=rules.plink_subject_filter.output.fam,
        id_map=rules.list_of_subjects_to_analyze.output.sample2subject,
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
    filename = checkpoints.selected_subjects.get(**wildcards).output[0]
    min_num_subjects = cfg.config.workflow_params.minimum_pop_subjects
    return (
        subject_qc_table.read(filename)
        .groupby("Ancestry")
        .size()
        .pipe(lambda x: x[x >= min_num_subjects])  # defaults to 50
        .index.tolist()
    )


# Split Subjects by Populations
rule split_subjects_by_population:
    input:
        rules.selected_subjects.output[0],
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
use rule keep_ids from plink as plink_split_subjects_by_population with:
    input:
        bed=rules.plink_convert_sample_id_to_subject_id.output.bed,
        bim=rules.plink_convert_sample_id_to_subject_id.output.bim,
        fam=rules.plink_convert_sample_id_to_subject_id.output.fam,
        to_keep="subject_level/{population}/subjects.txt",
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


use rule maf_filter from plink as plink_pop_maf_filter with:
    input:
        bed="subject_level/{population}/{prefix}.bed",
        bim="subject_level/{population}/{prefix}.bim",
        fam="subject_level/{population}/{prefix}.fam",
    params:
        maf="{maf}",
        out_prefix="subject_level/{population}/{prefix}_maf{maf}",
    output:
        bed=temp("subject_level/{population}/{prefix}_maf{maf}.bed"),
        bim=temp("subject_level/{population}/{prefix}_maf{maf}.bim"),
        fam=temp("subject_level/{population}/{prefix}_maf{maf}.fam"),
        nosex=temp("subject_level/{population}/{prefix}_maf{maf}.nosex"),
    log:
        "subject_level/{population}/{prefix}_maf{maf}.log",
    group:
        "{population}"


use rule ld from plink as plink_pop_ld with:
    input:
        bed="subject_level/{population}/{prefix}.bed",
        bim="subject_level/{population}/{prefix}.bim",
        fam="subject_level/{population}/{prefix}.fam",
    params:
        r2="{ld}",  # r2 threshold: currently 0.1
        out_prefix="subject_level/{population}/{prefix}_ld{ld}",
    output:
        to_keep=temp("subject_level/{population}/{prefix}_ld{ld}.prune.in"),  # Markers in approx. linkage equilibrium
        to_remove=temp("subject_level/{population}/{prefix}_ld{ld}.prune.out"),  # Markers in LD
        nosex=temp("subject_level/{population}/{prefix}_ld{ld}.nosex"),  # Markers in LD
    log:
        "subject_level/{population}/{prefix}_ld{ld}.log",
    group:
        "{population}"


use rule ld_filter from plink as plink_pop_ld_filter with:
    input:
        bed="subject_level/{population}/{prefix}.bed",
        bim="subject_level/{population}/{prefix}.bim",
        fam="subject_level/{population}/{prefix}.fam",
        to_keep="subject_level/{population}/{prefix}_ld{ld}.prune.in",
    params:
        out_prefix="subject_level/{population}/{prefix}_ld{ld}_pruned",
    output:
        bed=temp("subject_level/{population}/{prefix}_ld{ld}_pruned.bed"),
        bim=temp("subject_level/{population}/{prefix}_ld{ld}_pruned.bim"),
        fam=temp("subject_level/{population}/{prefix}_ld{ld}_pruned.fam"),
        nosex=temp("subject_level/{population}/{prefix}_ld{ld}_pruned.nosex"),
    log:
        "subject_level/{population}/{prefix}_ld{ld}_pruned.log",
    group:
        "{population}"


use rule keep_bfile from plink as plink_pop_keep_bfile with:
    input:
        bed="subject_level/{population}/{prefix}.bed",
        bim="subject_level/{population}/{prefix}.bim",
        fam="subject_level/{population}/{prefix}.fam",
    params:
        out_prefix="subject_level/{population}/{prefix}.keep",
    output:
        bed="subject_level/{population}/{prefix}.keep.bed",
        bim="subject_level/{population}/{prefix}.keep.bim",
        fam="subject_level/{population}/{prefix}.keep.fam",
        nosex="subject_level/{population}/{prefix}.keep.nosex",
    log:
        "subject_level/{population}/{prefix}.keep.log",
    group:
        "{population}"


use rule genome from plink as plink_pop_genome with:
    input:
        bed="subject_level/{population}/{prefix}.bed",
        bim="subject_level/{population}/{prefix}.bim",
        fam="subject_level/{population}/{prefix}.fam",
    params:
        ibd_min=cfg.config.software_params.ibd_pi_hat_min,
        ibd_max=cfg.config.software_params.ibd_pi_hat_max,
        out_prefix="subject_level/{population}/{prefix}",
    output:
        "subject_level/{population}/{prefix}.genome",
    group:
        "{population}"


rule subject_concordance_plink:
    input:
        "subject_level/{population}/{prefix}.genome",
    params:
        concordance_threshold=cfg.config.software_params.dup_concordance_cutoff,
        pi_hat_threshold=cfg.config.software_params.pi_hat_threshold,
    output:
        "subject_level/{population}/{prefix}.concordance.csv",
    group:
        "{population}"
    script:
        "../scripts/concordance_table.py"


def _population_concordance_files(wildcards):
    populations = _get_populations(wildcards)

    if not populations:
        return []

    return expand(
        "subject_level/{population}/subjects_maf{maf}_ld{ld}_pruned.keep.concordance.csv",
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
rule related_subjects:
    input:
        "subject_level/{{population}}/subjects_maf{maf}_ld{ld}_pruned.keep.concordance.csv".format(
            maf=cfg.config.software_params.maf_for_ibd,
            ld=cfg.config.software_params.ld_prune_r2,
        ),
    output:
        relatives="subject_level/{population}/relatives.csv",
        to_remove="subject_level/{population}/related_subjects_to_remove.txt",
    group:
        "{population}"
    script:
        "../scripts/related_subjects.py"


use rule remove_ids from plink as plink_remove_related_subjects with:
    input:
        bed="subject_level/{population}/subjects.bed",
        bim="subject_level/{population}/subjects.bim",
        fam="subject_level/{population}/subjects.fam",
        to_remove=rules.related_subjects.output.to_remove,
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
use rule bed_to_ped from plink as plink_convert_bed_to_ped with:
    input:
        bed="subject_level/{population}/{prefix}.bed",
        bim="subject_level/{population}/{prefix}.bim",
        fam="subject_level/{population}/{prefix}.fam",
    output:
        ped=temp("subject_level/{population}/{prefix}.ped"),
        map_=temp("subject_level/{population}/{prefix}.map"),
    params:
        out_prefix="subject_level/{population}/{prefix}",
    log:
        "subject_level/{population}/{prefix}.log",
    group:
        "{population}"


use rule convert from eigensoft as eigensoft_convert with:
    input:
        ped="subject_level/{population}/{prefix}.ped",
        map_="subject_level/{population}/{prefix}.map",
    output:
        par=temp("subject_level/{population}/{prefix}.convert.par"),
        gen=temp("subject_level/{population}/{prefix}.gen"),
        snp=temp("subject_level/{population}/{prefix}.snp"),
        ind=temp("subject_level/{population}/{prefix}.ind"),
    group:
        "{population}"


use rule smartpca from eigensoft as eigensoft_smartpca with:
    input:
        gen="subject_level/{population}/{prefix}.gen",
        snp="subject_level/{population}/{prefix}.snp",
        ind="subject_level/{population}/{prefix}.ind",
    output:
        par=temp("subject_level/{population}/{prefix}.pca.par"),
        eigenvec="subject_level/{population}/{prefix}.eigenvec",
    group:
        "{population}"


rule plot_pca:
    input:
        qc_table=rules.selected_subjects.output[0],
        eigenvec="subject_level/{{population}}/subjects_unrelated_maf{maf}_ld{ld}_pruned.keep.eigenvec".format(
            maf=cfg.config.software_params.maf_for_ibd,
            ld=cfg.config.software_params.ld_prune_r2,
        ),
    params:
        population="{population}",
    output:
        "subject_level/pca_plots/{population}.png",
    group:
        "{population}"
    script:
        "../scripts/plot_pca.py"


# Autosomal Heterozygosity (see plink_stats.smk 'plink_stats_het')
use rule het from plink as plink_controls_het with:
    input:
        bed="subject_level/{population}/{prefix}.bed",
        bim="subject_level/{population}/{prefix}.bim",
        fam="subject_level/{population}/{prefix}.fam",
    params:
        out_prefix="subject_level/{population}/{prefix}",
    output:
        "subject_level/{population}/{prefix}.het",
    group:
        "{population}"


rule plot_autosomal_heterozygosity:
    input:
        qc_table=rules.selected_subjects.output[0],
        het="subject_level/{population}/subjects.het",
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
        relatives="subject_level/{population}/relatives.csv",
        pca=rules.plot_pca.input.eigenvec,
        autosomal_het=rules.plot_autosomal_heterozygosity.input.het,
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
        subject_qc_table=rules.selected_subjects.output[0],
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
rule controls_per_population:
    input:
        rules.selected_subjects.output[0],
    output:
        "subject_level/{population}/control_list.txt",
    group:
        "{population}"
    run:
        (
            subject_qc_table.read(input[0])
            .query("Ancestry == @wildcards.population & case_control == 'Control'")
            .reindex(["Group_By_Subject_ID", "Group_By_Subject_ID"], axis=1)
            .to_csv(output[0], sep=" ", index=False, header=False)
        )


use rule keep_ids from plink as plink_split_controls_by_population with:
    input:
        bed="subject_level/{population}/subjects_unrelated.bed",
        bim="subject_level/{population}/subjects_unrelated.bim",
        fam="subject_level/{population}/subjects_unrelated.fam",
        to_keep=rules.controls_per_population.output[0],
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


use rule snps_only_filter from plink as plink_* with:
    input:
        bed="subject_level/{population}/{prefix}.bed",
        bim="subject_level/{population}/{prefix}.bim",
        fam="subject_level/{population}/{prefix}.fam",
    params:
        out_prefix="subject_level/{population}/{prefix}_snps",
    output:
        bed=temp("subject_level/{population}/{prefix}_snps.bed"),
        bim=temp("subject_level/{population}/{prefix}_snps.bim"),
        fam=temp("subject_level/{population}/{prefix}_snps.fam"),
        nosex=temp("subject_level/{population}/{prefix}_snps.nosex"),
    log:
        "subject_level/{population}/{prefix}_snps.log",
    group:
        "{population}"


use rule autosome_only_filter from plink as plink_* with:
    input:
        bed="subject_level/{population}/{prefix}.bed",
        bim="subject_level/{population}/{prefix}.bim",
        fam="subject_level/{population}/{prefix}.fam",
    params:
        out_prefix="subject_level/{population}/{prefix}_autosome",
    output:
        bed=temp("subject_level/{population}/{prefix}_autosome.bed"),
        bim=temp("subject_level/{population}/{prefix}_autosome.bim"),
        fam=temp("subject_level/{population}/{prefix}_autosome.fam"),
        nosex=temp("subject_level/{population}/{prefix}_autosome.nosex"),
    log:
        "subject_level/{population}/{prefix}_autosome.log",
    group:
        "{population}"


use rule hwe from plink as plink_* with:
    input:
        bed="subject_level/{population}/{prefix}.bed",
        bim="subject_level/{population}/{prefix}.bim",
        fam="subject_level/{population}/{prefix}.fam",
    params:
        out_prefix="subject_level/{population}/{prefix}",
    output:
        "subject_level/{population}/{prefix}.hwe",
    group:
        "{population}"


rule plot_hwe:
    input:
        "subject_level/{{population}}/controls_unrelated_maf{maf}_snps_autosome.keep.hwe".format(
            maf=cfg.config.software_params.maf_for_hwe,
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
    filename = checkpoints.selected_subjects.get(**wildcards).output[0]
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
