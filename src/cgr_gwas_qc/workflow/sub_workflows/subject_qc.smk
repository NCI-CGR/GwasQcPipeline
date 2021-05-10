from more_itertools import flatten

from cgr_gwas_qc import load_config
from cgr_gwas_qc.workflow.scripts import subject_qc_table


cfg = load_config()


################################################################################
# All Targets
################################################################################
rule all_subject_qc:
    input:
        "subject_level/population_qc.csv",
        "subject_level/concordance.csv",
        "subject_level/.population_plots.done",
        "subject_level/.control_plots.done",


################################################################################
# Imports
################################################################################
include: cfg.modules("common")
include: cfg.modules("plink_filters")
include: cfg.modules("plink_stats")


subworkflow sample_qc:
    snakefile: cfg.subworkflow("sample_qc")
    workdir: cfg.root.as_posix()


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
    script:
        "../scripts/subject_qc_table.py"


checkpoint selected_subjects:
    """Pull out subjects that have no analytic exclusions."""
    input:
        rules.subject_qc_table.output[0],
    output:
        "subject_level/selected_subjects.csv",
    run:
        (
            subject_qc_table.read(input[0])
            .query("not subject_analytic_exclusion")
            .to_csv(output[0], index=False)
        )


# Convert Sample Level data to Subject Level
rule files_for_plink_subject_filters:
    """Create mapping file for plink to convert Sample_IDs to Subject_IDs."""
    input:
        rules.selected_subjects.output[0],
    output:
        selected=temp("subject_level/selected_subjects.txt"),
        sample2subject=temp("subject_level/sample_to_subject.txt"),
    run:
        qc = subject_qc_table.read(input[0])
        qc.reindex(["Sample_ID", "Sample_ID"], axis=1).to_csv(
            output.selected, sep=" ", index=False, header=False
        )
        qc.reindex(
            ["Sample_ID", "Sample_ID", "Group_By_Subject_ID", "Group_By_Subject_ID"], axis=1,
        ).to_csv(output.sample2subject, sep=" ", index=False, header=False)


rule pull_subject_representative:
    input:
        bed=sample_qc("sample_level/call_rate_2/samples.bed"),
        bim=sample_qc("sample_level/call_rate_2/samples.bim"),
        fam=sample_qc("sample_level/call_rate_2/samples.fam"),
        selected=rules.files_for_plink_subject_filters.output.selected,
    params:
        out_prefix="subject_level/samples",
    output:
        bed="subject_level/samples.bed",
        bim="subject_level/samples.bim",
        fam="subject_level/samples.fam",
        nosex="subject_level/samples.nosex",
    log:
        "subject_level/samples.log",
    conda:
        cfg.conda("plink2")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--keep {input.selected} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule convert_Sample_ID_to_Subject_ID:
    input:
        bed=rules.pull_subject_representative.output.bed,
        bim=rules.pull_subject_representative.output.bim,
        fam=rules.pull_subject_representative.output.fam,
        sample2subject=rules.files_for_plink_subject_filters.output.sample2subject,
    params:
        out_prefix="subject_level/subjects",
    output:
        bed="subject_level/subjects.bed",
        bim="subject_level/subjects.bim",
        fam="subject_level/subjects.fam",
        nosex="subject_level/subjects.nosex",
    log:
        "subject_level/subjects.log",
    conda:
        cfg.conda("plink2")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--update-ids {input.sample2subject} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


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
    run:
        (
            subject_qc_table.read(input[0])
            .query("Ancestry == @wildcards.population")
            .reindex(["Group_By_Subject_ID", "Group_By_Subject_ID"], axis=1)
            .to_csv(output[0], sep=" ", index=False, header=False)
        )


rule split_plink_file_by_population:
    input:
        bed="subject_level/subjects.bed",
        bim="subject_level/subjects.bim",
        fam="subject_level/subjects.fam",
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
    wildcard_constraints:
        population="\w+",
    conda:
        cfg.conda("plink2")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--keep {input.to_keep} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


# Within Population Subject Concordance
def _population_concordance_files(wildcards):
    populations = _get_populations(wildcards)

    if not populations:
        return []

    return expand(
        "subject_level/{population}/subjects_maf{maf}_ld{ld}_pruned.concordance.csv",
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
rule related_subjects:
    input:
        "{{prefix}}/subjects_maf{maf}_ld{ld}_pruned.concordance.csv".format(
            maf=cfg.config.software_params.maf_for_ibd, ld=cfg.config.software_params.ld_prune_r2,
        ),
    output:
        relatives="{prefix}/relatives.csv",
        to_remove="{prefix}/related_subjects_to_remove.txt",
    script:
        "../scripts/related_subjects.py"


rule remove_related_subjects:
    input:
        bed="{prefix}/subjects.bed",
        bim="{prefix}/subjects.bim",
        fam="{prefix}/subjects.fam",
        to_remove=rules.related_subjects.output.to_remove,
    params:
        out_prefix="{prefix}/subjects_unrelated",
    output:
        bed="{prefix}/subjects_unrelated.bed",
        bim="{prefix}/subjects_unrelated.bim",
        fam="{prefix}/subjects_unrelated.fam",
        nosex="{prefix}/subjects_unrelated.nosex",
    log:
        "{prefix}/subjects_unrelated.log",
    conda:
        cfg.conda("plink2")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--remove {input.to_remove} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


# PCA
rule convert_bed_to_ped:
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
    conda:
        cfg.conda("plink2")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--recode "
        "--keep-allele-order "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


def eigensoft_config_inputs(wildcards):
    tool = wildcards.tool
    population = wildcards.population
    prefix = wildcards.prefix
    if tool == "pca":
        return {
            "gen": f"subject_level/{population}/{prefix}.gen",
            "snp": f"subject_level/{population}/{prefix}.snp",
            "ind": f"subject_level/{population}/{prefix}.ind",
        }

    # Default is to convert PED/MAP to EIGENSTRAT
    return {
        "gen": f"subject_level/{population}/{prefix}.ped",
        "snp": f"subject_level/{population}/{prefix}.map",
        "ind": f"subject_level/{population}/{prefix}.ped",
    }


def eigensoft_config_params(wildcards):
    tool = wildcards.tool
    population = wildcards.population
    prefix = wildcards.prefix

    if tool == "pca":
        return f"""\
        genotypename: subject_level/{population}/{prefix}.gen
        snpname: subject_level/{population}/{prefix}.snp
        indivname: subject_level/{population}/{prefix}.ind
        evecoutname: subject_level/{population}/{prefix}.eigenvec
        fastmode: YES
        """

    # Default is to convert PED/MAP to EIGENSTRAT
    return f"""\
    genotypename: subject_level/{population}/{prefix}.ped
    snpname: subject_level/{population}/{prefix}.map
    indivname: subject_level/{population}/{prefix}.ped
    outputformat: EIGENSTRAT
    genooutfilename: subject_level/{population}/{prefix}.gen
    snpoutfilename: subject_level/{population}/{prefix}.snp
    indoutfilename: subject_level/{population}/{prefix}.ind
    familynames: NO
    """


rule eigensoft_create_par:
    input:
        unpack(eigensoft_config_inputs),
    params:
        eigensoft_config_params,
    output:
        temp("subject_level/{population}/{prefix}.{tool}.par"),
    wildcard_constraints:
        tool="convert|pca",
    run:
        Path(output[0]).write_text(dedent(params[0]))


rule eigensoft_convert:
    input:
        ped=rules.convert_bed_to_ped.output.ped,
        map_=rules.convert_bed_to_ped.output.map_,
        par="subject_level/{population}/{prefix}.convert.par",
    output:
        gen=temp("subject_level/{population}/{prefix}.gen"),
        snp=temp("subject_level/{population}/{prefix}.snp"),
        ind=temp("subject_level/{population}/{prefix}.ind"),
    conda:
        cfg.conda("eigensoft")
    shell:
        "convertf -p {input.par}"


rule eigensoft_smartpca:
    input:
        gen=rules.eigensoft_convert.output.gen,
        snp=rules.eigensoft_convert.output.snp,
        ind=rules.eigensoft_convert.output.ind,
        par="subject_level/{population}/{prefix}.pca.par",
    output:
        gen="subject_level/{population}/{prefix}.eigenvec",
    conda:
        cfg.conda("eigensoft")
    shell:
        "smartpca -p {input.par}"


rule plot_pca:
    input:
        qc_table=rules.selected_subjects.output[0],
        eigenvec="subject_level/{{population}}/subjects_unrelated_maf{maf}_ld{ld}_pruned.eigenvec".format(
            maf=cfg.config.software_params.maf_for_ibd, ld=cfg.config.software_params.ld_prune_r2,
        ),
    params:
        population="{population}",
    output:
        "subject_level/pca_plots/{population}.png",
    script:
        "../scripts/plot_pca.py"


# Autosomal Heterozygosity (see plink_stats.smk 'plink_stats_het')
rule plot_autosomal_heterozygosity:
    input:
        qc_table=rules.selected_subjects.output[0],
        het="subject_level/{population}/subjects.het",
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
        relatives="subject_level/{population}/relatives.csv",
        pca=rules.plot_pca.input.eigenvec,
        autosomal_het=rules.plot_autosomal_heterozygosity.input.het,
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
        subject_qc_table=rules.selected_subjects.output[0],
        population_qc_tables=_population_qc_tables,
    output:
        "subject_level/population_qc.csv",
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
    run:
        (
            subject_qc_table.read(input[0])
            .query("Ancestry == @wildcards.population & case_control == 'Control'")
            .reindex(["Group_By_Subject_ID", "Group_By_Subject_ID"], axis=1)
            .to_csv(output[0], sep=" ", index=False, header=False)
        )


rule plink_split_controls:
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
    log:
        "subject_level/{population}/controls_unrelated.log",
    conda:
        cfg.conda("plink2")
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
    shell:
        "plink "
        "--bed {input.bed} "
        "--bim {input.bim} "
        "--fam {input.fam} "
        "--keep {input.to_keep} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


rule plot_hwe:
    input:
        "subject_level/{{population}}/controls_unrelated_maf{maf}_snps_autosome_cleaned.hwe".format(
            maf=cfg.config.software_params.maf_for_hwe,
        ),
    params:
        population="{population}",
    output:
        "subject_level/hwe_plots/{population}.png",
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
