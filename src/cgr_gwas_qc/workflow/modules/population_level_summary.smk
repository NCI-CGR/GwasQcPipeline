from pathlib import Path

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
        "sample_level/sample_qc.csv",
    params:
        threshold=cfg.config.workflow_params.minimum_pop_subjects,
    output:
        directory("population_level/subject_lists"), #"population_level/subject_lists/{population}.txt"
    run:
        output_path = Path(output[0])
        output_path.mkdir(exist_ok=True, parents=True)

        df = pd.read_csv(input[0]).query("is_subject_representative")
        for pop_, grp in df.groupby("Ancestry"):
            if grp.shape[0] < params.threshold:
                # Too few subjects to analyze population
                continue

            (  # Save a list of subjects for each population
                grp.assign(Subject_ID2=lambda x: x.Group_By_Subject_ID)
                .reindex(["Group_By_Subject_ID", "Subject_ID2"], axis=1)
                .to_csv(output_path / f"{pop_}.txt", sep=" ", index=False, header=False)
            )

        if len(list(output_path.glob("*.txt"))) == 0:
            (output_path / "no_populations.txt").touch()


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


rule population_qc_table:
    input:
        relatives="population_level/{{population}}/subjects_relatives_pi_hat_gt{pi}.csv".format(
            pi=cfg.config.software_params.pi_hat_threshold
        ),
        pca="population_level/{{population}}/subjects_unrelated{pi}_maf{maf}_ld{ld}_pruned.eigenvec".format(
            pi=cfg.config.software_params.pi_hat_threshold,
            maf=cfg.config.software_params.maf_for_ibd,
            ld=cfg.config.software_params.ld_prune_r2,
        ),
        autosomal_het="population_level/{population}/subjects.het",
    params:
        population="{population}",
    output:
        "population_level/{population}/qc.csv",
    script:
        "../scripts/population_qc_table.py"


def required_population_results(wildcards):
    """Decide what populations to analyze.

    If a population has fewer than `workflow_params.minimum_pop_subjects`
    subjects than ignore.
    """
    _ = checkpoints.subjects_per_population.get(**wildcards).output[0]
    populations = glob_wildcards("population_level/subject_lists/{population}.txt").population

    if populations == ["no_populations"]:
        return []

    return expand("population_level/{population}/qc.csv", population=populations)


rule per_population_qc_done:
    input:
        required_population_results,
    output:
        "population_level/per_population_qc.done",
    shell:
        "echo {input} | xargs printf '%s\n' > {output[0]}"


################################################################################
# Population Level Analysis (Controls Only)
#   - HWE
################################################################################
checkpoint controls_per_population:
    input:
        "sample_level/sample_qc.csv",
    params:
        threshold=cfg.config.workflow_params.control_hwp_threshold,
    output:
        directory("population_level/controls_lists"),
    run:
        output_path = Path(output[0])
        output_path.mkdir(exist_ok=True, parents=True)

        df = pd.read_csv(input[0]).query("is_subject_representative & case_control == 'control'")
        for pop_, grp in df.groupby("Ancestry"):
            if grp.shape[0] < params.threshold:
                # Too few controls to analyze population
                continue

            (  # Save a list of subjects for each population
                grp.assign(Subject_ID2=lambda x: x.Group_By_Subject_ID)
                .reindex(["Group_By_Subject_ID", "Subject_ID2"], axis=1)
                .to_csv(output_path / f"{pop_}.txt", sep=" ", index=False, header=False)
            )

        if len(list(output_path.glob("*.txt"))) == 0:
            (output_path / "no_controls.txt").touch()


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

    if populations == ["no_controls"]:
        return []

    maf = cfg.config.software_params.maf_for_hwe
    pi = cfg.config.software_params.pi_hat_threshold

    return expand(  # HWE
        "population_level/{population}/controls_unrelated{pi}_maf{maf}_snps_autosome_cleaned.hwe",
        population=populations,
        maf=maf,
        pi=pi,
    )


rule per_population_controls_qc_done:
    input:
        required_population_controls,
    output:
        "population_level/per_population_controls_qc.done",
    shell:
        "echo {input} | xargs printf '%s\n' > {output[0]}"


################################################################################
# Summary of Population and Controls
################################################################################
rule agg_population_qc_table:
    input:
        sample_qc_table="sample_level/sample_qc.csv",
        population_qc_tables=rules.per_population_qc_done.output[0],
    output:
        "population_level/population_qc.csv",
    script:
        "../scripts/agg_population_qc_tables.py"
