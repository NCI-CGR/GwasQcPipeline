from pathlib import Path

from more_itertools import flatten

from cgr_gwas_qc import load_config
from cgr_gwas_qc.workflow.scripts import subject_qc_table
import shutil

cfg = load_config()


localrules:
    all_subject_qc,
    subject_qc_table,
    selected_Subject_IDs,
    population_checkpoint,
    population_Subject_IDs,
    population_level_concordance_plink,
    agg_population_concordance,
    per_population_qc_table,
    agg_population_qc_tables,
    plot_pca,
    plot_autosomal_heterozygosity,
    agg_population_plots,
    population_controls_checkpoint,
    population_controls_Subject_IDs,
    plot_hwe,
    agg_control_plots,    

wildcard_constraints:
    population="[\w_]+",


################################################################################
# Subject QC Targets
################################################################################
targets = [
    "subject_level/subject_qc.csv",
    "subject_level/samples.bed",
    "subject_level/samples.bim",
    "subject_level/samples.fam",
    "subject_level/subjects.bed",
    "subject_level/subjects.bim",
    "subject_level/subjects.fam",
    "subject_level/concordance.csv",
    "subject_level/population_qc.csv",
    "subject_level/.population_plots.done",
    "subject_level/.control_plots.done",
    "delivery/gwas.assoc",
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
    output:
        "subject_level/subject_qc.csv",
    script:
        "../scripts/subject_qc_table.py"


# Convert Sample Level data to Subject Level
rule selected_Subject_IDs:
    """Create mapping file for plink to convert Sample_IDs to Subject_IDs."""
    input:
        rules.subject_qc_table.output[0],
    output:
        selected=temp("subject_level/selected_Subject_IDs.txt"),
        sample2subject=temp("subject_level/sample_to_subject.txt"),
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


# -------------------------------------------------------------------------------
# Population Level Analysis
# -------------------------------------------------------------------------------
checkpoint population_checkpoint:
    input:
        rules.subject_qc_table.output[0],
    params:
        min_num_subjects=cfg.config.workflow_params.minimum_pop_subjects,
    output:
        directory("subject_level/populations"),
    run:
        populations = (
            subject_qc_table.read(input[0])
            .groupby("Ancestry")
            .size()
            .pipe(lambda x: x[x >= params.min_num_subjects])
            .index.tolist()
        )

        path = Path(output[0])
        path.mkdir(exist_ok=True, parents=True)

        for population in populations:
            (path / population).touch()


def _get_populations(wildcards):
    checkpoint_output = checkpoints.population_checkpoint.get(**wildcards).output[0]
    return [
        x
        for x in glob_wildcards(Path(checkpoint_output, "{population}")).population
        if not x.startswith(".snakemake")
    ]


# Split Subjects by Populations
rule population_Subject_IDs:
    input:
        rules.subject_qc_table.output[0],
    output:
        "subject_level/{population}/subjects.txt",
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


use rule ld from plink as population_level_ld_estimate with:
    input:
        bed=rules.population_level_maf_filter.output.bed,
        bim=rules.population_level_maf_filter.output.bim,
        fam=rules.population_level_maf_filter.output.fam,
    params:
        r2="{ld}",  # r2 threshold: currently 0.1
        out_prefix="subject_level/{population}/subjects_maf{maf}_ld{ld}_estimate",
    output:
        # Markers in approx. linkage equilibrium
        to_keep=temp("subject_level/{population}/subjects_maf{maf}_ld{ld}_estimate.prune.in"),
        # Markers in LD
        to_remove=temp("subject_level/{population}/subjects_maf{maf}_ld{ld}_estimate.prune.out"),
        nosex=temp("subject_level/{population}/subjects_maf{maf}_ld{ld}_estimate.nosex"),
    log:
        "subject_level/{population}/subjects_maf{maf}_ld{ld}_estimate.log",


use rule ld_filter from plink as population_level_ld_pruned with:
    input:
        bed=rules.population_level_maf_filter.output.bed,
        bim=rules.population_level_maf_filter.output.bim,
        fam=rules.population_level_maf_filter.output.fam,
        to_keep=rules.population_level_ld_estimate.output.to_keep,
    params:
        out_prefix="subject_level/{population}/subjects_maf{maf}_ld{ld}",
    output:
        bed="subject_level/{population}/subjects_maf{maf}_ld{ld}.bed",
        bim="subject_level/{population}/subjects_maf{maf}_ld{ld}.bim",
        fam="subject_level/{population}/subjects_maf{maf}_ld{ld}.fam",
        nosex="subject_level/{population}/subjects_maf{maf}_ld{ld}.nosex",
    log:
        "subject_level/{population}/subjects_maf{maf}_ld{ld}.log",


use rule genome from plink as population_level_ibd with:
    input:
        bed=rules.population_level_ld_pruned.output.bed,
        bim=rules.population_level_ld_pruned.output.bim,
        fam=rules.population_level_ld_pruned.output.fam,
    params:
        ibd_min=cfg.config.software_params.ibd_pi_hat_min,
        ibd_max=cfg.config.software_params.ibd_pi_hat_max,
        out_prefix="subject_level/{population}/subjects_maf{maf}_ld{ld}_ibd",
    output:
        "subject_level/{population}/subjects_maf{maf}_ld{ld}_ibd.genome",


rule population_level_concordance_plink:
    input:
        rules.population_level_ibd.output[0],
    params:
        concordance_threshold=cfg.config.software_params.dup_concordance_cutoff,
        pi_hat_threshold=cfg.config.software_params.pi_hat_threshold,
    output:
        "subject_level/{population}/subjects_maf{maf}_ld{ld}.concordance.csv",
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


use rule ld from plink as population_level_unrelated_ld_estimate with:
    input:
        bed=rules.population_level_unrelated_maf_filter.output.bed,
        bim=rules.population_level_unrelated_maf_filter.output.bim,
        fam=rules.population_level_unrelated_maf_filter.output.fam,
    params:
        r2="{ld}",  # r2 threshold: currently 0.1
        out_prefix="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_estimate",
    output:
        # Markers in approx. linkage equilibrium
        to_keep=temp(
            "subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_estimate.prune.in"
        ),
        # Markers in LD
        to_remove=temp(
            "subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_estimate.prune.out"
        ),
        nosex=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_estimate.nosex"),
    log:
        "subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_estimate.log",


use rule ld_filter from plink as population_level_unrelated_ld_pruned with:
    input:
        bed=rules.population_level_unrelated_maf_filter.output.bed,
        bim=rules.population_level_unrelated_maf_filter.output.bim,
        fam=rules.population_level_unrelated_maf_filter.output.fam,
        to_keep=rules.population_level_unrelated_ld_estimate.output.to_keep,
    params:
        out_prefix="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}",
    output:
        bed="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.bed",
        bim="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.bim",
        fam="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.fam",
        nosex="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.nosex",
    log:
        "subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.log",


use rule bed_to_ped from plink as population_level_unrelated_bed_to_ped with:
    input:
        bed=rules.population_level_unrelated_ld_pruned.output.bed,
        bim=rules.population_level_unrelated_ld_pruned.output.bim,
        fam=rules.population_level_unrelated_ld_pruned.output.fam,
    output:
        ped=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_ped_to_bed.ped"),
        map_=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_ped_to_bed.map"),
    params:
        out_prefix="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_ped_to_bed",
    log:
        "subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_ped_to_bed.log",


rule trim_ids:
    """EIGENSOFT convert requires sample/snp IDs are <39 characters."""
    input:
        ped=rules.population_level_unrelated_bed_to_ped.output.ped,
        map_=rules.population_level_unrelated_bed_to_ped.output.map_,
    output:
        ped=temp(
            "subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_ped_to_bed_trimmed.ped"
        ),
        map_=temp(
            "subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}_ped_to_bed_trimmed.map"
        ),
    script:
        "../scripts/trim_ped_map_ids.py"


use rule convert from eigensoft as population_level_unrelated_convert_to_eigensoft with:
    input:
        ped=rules.trim_ids.output.ped,
        map_=rules.trim_ids.output.map_,
    output:
        par=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.convert.par"),
        gen=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.gen"),
        snp=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.snp"),
        ind=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.ind"),


use rule smartpca from eigensoft as population_level_unrelated_smartpca with:
    input:
        gen=rules.population_level_unrelated_convert_to_eigensoft.output.gen,
        snp=rules.population_level_unrelated_convert_to_eigensoft.output.snp,
        ind=rules.population_level_unrelated_convert_to_eigensoft.output.ind,
    output:
        par=temp("subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.pca.par"),
        eigenvec="subject_level/{population}/subjects_unrelated_maf{maf}_ld{ld}.eigenvec",


rule plot_pca:
    input:
        qc_table=rules.subject_qc_table.output[0],
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


rule plot_autosomal_heterozygosity:
    input:
        qc_table=rules.subject_qc_table.output[0],
        het=rules.population_level_autosomal_heterozygosity.output[0],
    params:
        population="{population}",
        threshold=cfg.config.software_params.autosomal_het_threshold,
    output:
        "subject_level/autosomal_heterozygosity_plots/{population}.png",
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
    script:
        "../scripts/population_qc_table.py"


def _population_qc_tables(wildcards):
    populations = _get_populations(wildcards)

    if not populations:
        return []

    return expand(rules.per_population_qc_table.output[0], population=populations)


rule agg_population_qc_tables:
    input:
        subject_qc_table=rules.subject_qc_table.output[0],
        population_qc_tables=_population_qc_tables,
    output:
        "subject_level/population_qc.csv",
    script:
        "../scripts/agg_population_qc_tables.py"

# gwas  (see plink_stats.smk 'gwas')
use rule gwas from plink as gwas with:
    input:
        bed="sample_level/samples.bed",
        bim="sample_level/samples.bim",
        fam="sample_level/samples.fam",
    params:
        out_prefix="delivery/gwas",
    output:
        "delivery/gwas.assoc",

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
checkpoint population_controls_checkpoint:
    input:
        rules.subject_qc_table.output[0],
    params:
        min_num_subjects=cfg.config.workflow_params.control_hwp_threshold,
    output:
        directory("subject_level/controls"),
    run:
        # issue 221  
        populations1 = (
            subject_qc_table.read(input[0])
            .query("case_control == 'Control'")
            .groupby("Ancestry")
            .size()
            .pipe(lambda x: x[x >= params.min_num_subjects])
            .index.tolist()
        )
        print(populations1)
        populations2 = (
            subject_qc_table.read(input[0])
            .query("case_control == 'Case'")
            .groupby("Ancestry")
            .size()
            .pipe(lambda x: x[x >= params.min_num_subjects])
            .index.tolist()
        )
        print(populations2)
        populations3 = (
            subject_qc_table.read(input[0])
            .query("case_control == 'Unknown'")
            .groupby("Ancestry")
            .size()
            .pipe(lambda x: x[x >= params.min_num_subjects])
            .index.tolist()
        )
        print(populations3)
        populations = populations1 + populations2
        populations = list(set(populations))
        print(populations)

        path = Path(output[0])
        path.mkdir(exist_ok=True, parents=True)

        for population in populations:
            (path / population).touch()


def _get_controls(wildcards):
    checkpoint_output = checkpoints.population_controls_checkpoint.get(**wildcards).output[0]
    return [
        x
        for x in glob_wildcards(Path(checkpoint_output, "{population}")).population
        if not x.startswith(".snakemake")
    ]


rule population_controls_Subject_IDs:
    input:
        rules.subject_qc_table.output[0],
    output:
        "subject_level/{population}/controls.txt",
    params:
        "subject_level/{population}/test_controls.txt",
        "subject_level/{population}/test_case.txt",
        "subject_level/{population}/test_controlscase.txt",
        "subject_level/{population}/unknown.txt",
        cfg.config.workflow_params.control_hwp_threshold,
        
    run:
        # issue 221,235 fix
        (
            subject_qc_table.read(input[0])
            .query("Ancestry == @wildcards.population & case_control == 'Control'")
            .reindex(["Group_By_Subject_ID", "Group_By_Subject_ID"], axis=1)
            .to_csv(params[0], sep=" ", index=False, header=False)
        )
        (
            subject_qc_table.read(input[0])
            .query("Ancestry == @wildcards.population & case_control == 'Case'")
            .reindex(["Group_By_Subject_ID", "Group_By_Subject_ID"], axis=1)
            .to_csv(params[1], sep=" ", index=False, header=False)
        )
        (
           subject_qc_table.read(input[0])
           .query("Ancestry == @wildcards.population & (case_control == 'Case'| case_control == 'Control')")
           .reindex(["Group_By_Subject_ID", "Group_By_Subject_ID"], axis=1)
           .to_csv(params[2], sep=" ", index=False, header=False)
        )
        (
           subject_qc_table.read(input[0])
           .query("Ancestry == @wildcards.population & case_control == 'Unknown'")
           .reindex(["Group_By_Subject_ID", "Group_By_Subject_ID"], axis=1)
           .to_csv(params[3], sep=" ", index=False, header=False)
        )

        min_num_subjects = int(params[4])
        print("minimum min_num_subjects:", min_num_subjects)
        # Controls 
        print(params[0])
        num_lines0 = sum(1 for line in open(params[0]))
        print(num_lines0)
        # Case
        print(params[1])
        num_lines1 = sum(1 for line in open(params[1]))
        print(num_lines1)
        # Controls + Case
        print(params[2])
        num_lines2 = sum(1 for line in open(params[2]))
        print(num_lines2)
        # Unknown
        print(params[3])
        num_lines3 = sum(1 for line in open(params[3]))
        print(num_lines3)

      
        if num_lines0 <  min_num_subjects:
            print(num_lines0, "controls samples less than", min_num_subjects, "threshold, using control and case samples for HWE")
            shutil.copyfile(params[2], output[0])
        else:
            print(num_lines0, "control samples present and above the ", min_num_subjects, "threshold for HWE")
            shutil.copyfile(params[0], output[0])

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
        "{population}_controls_filter"


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
        "{population}_controls_filter"


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
        "{population}_controls_filter"


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
        "{population}_controls_filter"


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
    script:
        "../scripts/plot_hwe.py"


def _control_plots(wildcards):
    populations = _get_controls(wildcards)

    if not populations:
        return []

    return expand(rules.plot_hwe.output[0], population=populations)


rule agg_control_plots:
    input:
        _control_plots,
    output:
        touch("subject_level/.control_plots.done"),
