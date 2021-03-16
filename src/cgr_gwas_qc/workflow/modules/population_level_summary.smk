import pandas as pd
from more_itertools import flatten


include: cfg.modules("common.smk")


################################################################################
# Population Level Analysis
#   - PCA
#   - Autosomal Heterozygosity
#   - IBS/IBD
################################################################################
checkpoint subjects_per_population:
    input:
        "sample_level/qc_summary.csv",
    params:
        threshold=cfg.config.workflow_params.minimum_pop_subjects,
    output:
        directory("population_level/subject_lists"), #"population_level/subject_lists/{population}.txt"
    run:
        df = pd.read_csv(input[0]).query("Subject_Representative")
        for pop_, grp in df.groupby("Ancestry"):
            if grp.shape[0] < params.threshold:
                # Too few subjects to analyze population
                continue

            pop_path = Path(output[0])
            pop_path.mkdir(exist_ok=True, parents=True)

            (  # Save a list of subjects for each population
                grp.assign(Subject_ID2=lambda x: x.Group_By_Subject_ID)
                .reindex(["Group_By_Subject_ID", "Subject_ID2"], axis=1)
                .to_csv(pop_path / f"{pop_}.txt", sep=" ", index=False, header=False)
            )


rule plink_split_population:
    input:
        bed="subject_level/subjects.bed",
        bim="subject_level/subjects.bim",
        fam="subject_level/subjects.fam",
        to_keep="population_level/subject_lists/{population}.txt",
    params:
        out_prefix="population_level/{population}/subjects",
    output:
        bed="population_level/{population}/subjects.bed",
        bim="population_level/{population}/subjects.bim",
        fam="population_level/{population}/subjects.fam",
        nosex="population_level/{population}/subjects.nosex",
    log:
        "population_level/{population}/subjects.log",
    wildcard_constraints:
        population="\w+",
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
        "--keep {input.to_keep} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


def required_population_results(wildcards):
    """Decide what populations to analyze.

    If a population has fewer than `workflow_params.minimum_pop_subjects`
    subjects than ignore.
    """
    _ = checkpoints.subjects_per_population.get(**wildcards).output[0]
    populations = glob_wildcards("population_level/subject_lists/{population}.txt").population

    pi = cfg.config.software_params.pi_hat_threshold
    maf = cfg.config.software_params.maf_for_ibd
    ld = cfg.config.software_params.ld_prune_r2

    return flatten(
        [
            expand(  # PCA
                "population_level/{population}/subjects_unrelated{pi}_maf{maf}_ld{ld}_pruned.eigenvec",
                population=populations,
                pi=pi,
                maf=maf,
                ld=ld,
            ),
            expand(  # IBS/IBD
                "population_level/{population}/subjects_maf{maf}_ld{ld}_pruned.genome",
                population=populations,
                maf=maf,
                ld=ld,
            ),
            expand(  # Autosomal Heterozygosity
                "population_level/{population}/subjects.het",
                population=populations,
                maf=maf,
                ld=ld,
            ),
        ]
    )


rule phony_population_results:
    input:
        required_population_results,
    output:
        "population_level/results.done",
    shell:
        "echo {input} | xargs printf '%s\n' > {output[0]}"


################################################################################
# Population Level Analysis (Controls Only)
#   - HWE
################################################################################
checkpoint controls_per_population:
    input:
        "sample_level/qc_summary.csv",
    params:
        threshold=cfg.config.workflow_params.control_hwp_threshold,
    output:
        directory("population_level/controls_lists"),
    run:
        df = pd.read_csv(input[0]).query("Subject_Representative & `Case/Control_Status` == 0")
        for pop_, grp in df.groupby("Ancestry"):
            if grp.shape[0] < params.threshold:
                # Too few controls to analyze population
                continue

            pop_path = Path(output[0])
            pop_path.mkdir(exist_ok=True, parents=True)

            (  # Save a list of subjects for each population
                grp.assign(Subject_ID2=lambda x: x.Group_By_Subject_ID)
                .reindex(["Group_By_Subject_ID", "Subject_ID2"], axis=1)
                .to_csv(pop_path / f"{pop_}.txt", sep=" ", index=False, header=False)
            )


rule plink_split_controls:
    input:
        bed="population_level/{population}/subjects_unrelated{pi}.bed",
        bim="population_level/{population}/subjects_unrelated{pi}.bim",
        fam="population_level/{population}/subjects_unrelated{pi}.fam",
        to_keep="population_level/controls_lists/{population}.txt",
    params:
        out_prefix="population_level/{population}/controls_unrelated{pi}",
    output:
        bed="population_level/{population}/controls_unrelated{pi}.bed",
        bim="population_level/{population}/controls_unrelated{pi}.bim",
        fam="population_level/{population}/controls_unrelated{pi}.fam",
        nosex="population_level/{population}/controls_unrelated{pi}.nosex",
    log:
        "population_level/{population}/controls_unrelated{pi}.log",
    wildcard_constraints:
        pi="[01].\d+",
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
        "--keep {input.to_keep} "
        "--make-bed "
        "--threads {threads} "
        "--memory {resources.mem_mb} "
        "--out {params.out_prefix}"


def required_population_controls(wildcards):
    """Decide what populations controls to analyze.

    If population controls have fewer than `workflow_params.control_hwp_threshold`
    subjects than ignore.
    """
    _ = checkpoints.controls_per_population.get(**wildcards).output[0]
    populations = glob_wildcards("population_level/controls_lists/{population}.txt").population

    maf = cfg.config.software_params.maf_for_hwe
    pi = cfg.config.software_params.pi_hat_threshold

    return expand(  # HWE
        "population_level/{population}/controls_unrelated{pi}_maf{maf}_snps_autosome_cleaned.hwe",
        population=populations,
        maf=maf,
        pi=pi,
    )


rule phony_population_controls:
    input:
        required_population_controls,
    output:
        "population_level/controls.done",
    shell:
        "echo {input} | xargs printf '%s\n' > {output[0]}"


################################################################################
# Summary of Population and Controls
################################################################################
rule population_qc_table:
    input:
        sample_qc="sample_level/qc_summary.csv",
        populations=rules.phony_population_results.output[0],
        controls=rules.phony_population_controls.output[0],
    output:
        "population_level/population_qc.csv",
    script:
        "../scripts/population_qc_table.py"
