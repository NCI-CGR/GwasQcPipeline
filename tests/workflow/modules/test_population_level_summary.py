import re
import shutil

import pytest

from cgr_gwas_qc import load_config
from cgr_gwas_qc.testing import chdir, file_hashes_equal, run_snakemake, sorted_file_equal
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_subjects_per_population(tmp_path, qc_summary):
    # GIVEN: real data config and qc summary table
    data_cache = (
        RealData(tmp_path)
        .copy_sample_sheet()
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("population_level_summary.smk")

            rule all:
                input:
                    "population_agg.txt",

            def _input(wildcards):
                _ = checkpoints.subjects_per_population.get(**wildcards).output[0]
                subject_list = "population_level/subject_lists/{population}.txt"
                return expand(
                    subject_list,
                    population=glob_wildcards(subject_list).population
                )

            rule agg:
                input:
                    _input
                output:
                    "population_agg.txt"
                shell:
                    "echo {input} > {output[0]}"

            """
        )
    )
    (tmp_path / "sample_level").mkdir()
    shutil.copyfile(qc_summary, tmp_path / "sample_level/qc_summary.csv")

    # WHEN: run snakemake to generate population keep lists
    run_snakemake(tmp_path)

    # THEN: The sorted keep lists should match legacy
    assert sorted_file_equal(
        tmp_path / "population_level/subject_lists/EUR.txt",
        data_cache / "production_outputs/split_by_pop/EUR.keep.txt",
    )


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_plink_split_population(tmp_path, conda_envs):
    # GIVEN: real data config, subject level data, and the European subject list
    # NOTE: we have to use EUR b/c other populations don't have enough subjects
    # in test data
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy_sample_sheet()
        .copy(
            "production_outputs/split_by_pop/EUR.keep.txt",
            "population_level/subject_lists/EUR.txt",
        )
        .copy("production_outputs/subject_level/subjects.bed", "subject_level/subjects.bed",)
        .copy("production_outputs/subject_level/subjects.bim", "subject_level/subjects.bim",)
        .copy("production_outputs/subject_level/subjects.fam", "subject_level/subjects.fam",)
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("population_level_summary.smk")

            rule all:
                input:
                    "population_level/EUR/subjects.bed",
                    "population_level/EUR/subjects.bim",
                    "population_level/EUR/subjects.fam",
            """
        )
    )
    # WHEN: run snakemake to generate European data subsets
    run_snakemake(tmp_path)

    # THEN: the plink data sets match legcay
    assert file_hashes_equal(
        tmp_path / "population_level/EUR/subjects.bed",
        data_cache / "production_outputs/split_by_pop/EUR_subjects.bed",
    )

    assert file_hashes_equal(
        tmp_path / "population_level/EUR/subjects.bim",
        data_cache / "production_outputs/split_by_pop/EUR_subjects.bim",
    )

    assert file_hashes_equal(
        tmp_path / "population_level/EUR/subjects.fam",
        data_cache / "production_outputs/split_by_pop/EUR_subjects.fam",
    )


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
@pytest.mark.slow
def test_phony_population_results(tmp_path, conda_envs, qc_summary):
    # GIVEN: real data config, qc_summary table, all of the inputs to generate
    # the summary table, and subject level plink data sets.
    # NOTE: I had to include the qc_summary inputs to get snakemake to run the
    # check point correctly.
    conda_envs.copy_env("plink2", tmp_path)
    conda_envs.copy_env("eigensoft", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy_sample_sheet()
        .copy("production_outputs/subject_level/subjects.bed", "subject_level/subjects.bed",)
        .copy("production_outputs/subject_level/subjects.bim", "subject_level/subjects.bim",)
        .copy("production_outputs/subject_level/subjects.fam", "subject_level/subjects.fam",)
        .make_config(workflow_params={"subject_id_to_use": "PI_Subject_ID"})
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("common.smk")
            include: cfg.modules("plink_stats.smk")
            include: cfg.modules("plink_filters.smk")
            include: cfg.modules("subject_level_summary.smk")
            include: cfg.modules("population_level_summary.smk")

            rule all:
                input:
                    "population_level/results.done",

            """
        )
    )
    (tmp_path / "sample_level").mkdir()
    shutil.copyfile(qc_summary, tmp_path / "sample_level/qc_summary.csv")

    # WHEN: run snakemake to get all population level and all population-control level results
    run_snakemake(tmp_path)

    # THEN:
    with chdir(tmp_path):
        cfg = load_config()
        maf = cfg.config.software_params.maf_for_ibd
        ld = cfg.config.software_params.ld_prune_r2
        pi = cfg.config.software_params.pi_hat_threshold

    assert not (tmp_path / "population_level/AFR").exists()  # Should not run b/c too few subjects
    assert not (tmp_path / "population_level/ASN").exists()  # Should not run b/c too few subjects

    assert file_hashes_equal(
        tmp_path / "population_level/EUR/subjects.bed",
        data_cache / "production_outputs/split_by_pop/EUR_subjects.bed",
    )

    # NOTE: The eigenvec comparison is flaky b/c signs change. I am scrubbing
    # that for the comparisons.
    def clean_eigenvec(pth):
        return re.sub("-", "", re.sub(r"\s+", " ", pth.read_text()))

    assert clean_eigenvec(
        tmp_path / f"population_level/EUR/subjects_unrelated{pi}_maf{maf}_ld{ld}_pruned.eigenvec"
    ) == clean_eigenvec(data_cache / "production_outputs/pca/EUR_subjects.eigenvec")

    assert file_hashes_equal(
        tmp_path / "population_level/EUR/subjects.het",
        data_cache / "production_outputs/autosomal_heterozygosity/EUR_subjects_qc.het",
    )


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_controls_per_population(tmp_path, qc_summary):
    # GIVEN: real data config and the qc_summary table
    # NOTE: we have to use EUR b/c other populations don't have enough subjects
    # in test data
    data_cache = (
        RealData(tmp_path)
        .copy_sample_sheet()
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("population_level_summary.smk")

            rule all:
                input:
                    "population_agg.txt",

            def _input(wildcards):
                _ = checkpoints.controls_per_population.get(**wildcards).output[0]
                controls_list = "population_level/controls_lists/{population}.txt"
                return expand(
                    controls_list,
                    population=glob_wildcards(controls_list).population
                )

            rule agg:
                input:
                    _input
                output:
                    "population_agg.txt"
                shell:
                    "echo {input} > {output[0]}"

            """
        )
    )
    (tmp_path / "sample_level").mkdir()
    shutil.copyfile(qc_summary, tmp_path / "sample_level/qc_summary.csv")

    # WHEN: run snakemake to generate population level control lists
    run_snakemake(tmp_path)

    # THEN: The lis of control subjects per population matches legacy
    assert sorted_file_equal(
        tmp_path / "population_level/controls_lists/EUR.txt",
        data_cache / "production_outputs/HWP/EUR_controls.txt",
    )


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_plink_split_controls(tmp_path, conda_envs):
    # GIVEN: real data config, a list of European control subjects and the European subjects
    # NOTE: we have to use EUR b/c other populations don't have enough subjects
    # in test data
    conda_envs.copy_env("plink2", tmp_path)

    pi = 0.2
    maf = 0.05

    data_cache = (
        RealData(tmp_path)
        .copy_sample_sheet()
        .copy("production_outputs/HWP/EUR_controls.txt", "population_level/controls_lists/EUR.txt",)
        .copy(
            "production_outputs/split_by_pop/EUR_subjects.bed",
            f"population_level/EUR/subjects_unrelated{pi}.bed",
        )
        .copy(
            "production_outputs/split_by_pop/EUR_subjects.bim",
            f"population_level/EUR/subjects_unrelated{pi}.bim",
        )
        .copy(
            "production_outputs/split_by_pop/EUR_subjects.fam",
            f"population_level/EUR/subjects_unrelated{pi}.fam",
        )
        .make_config(software_params={"maf_for_hwe": maf, "pi_hat_threshold": pi})
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("plink_filters.smk")
            include: cfg.modules("population_level_summary.smk")

            rule all:
                input:
                    expand(
                        "population_level/EUR/controls_unrelated{pi}_maf{maf}_snps_autosome_cleaned.{ext}",
                        maf=cfg.config.software_params.maf_for_hwe,
                        pi=cfg.config.software_params.pi_hat_threshold,
                        ext=["bed", "bim", "fam"]
                    )
            """
        )
    )

    # WHEN: run snakemake to generate the European control data sets.
    run_snakemake(tmp_path)

    # THEN: The European control data sets should match legacy
    assert file_hashes_equal(
        tmp_path
        / f"population_level/EUR/controls_unrelated{pi}_maf{maf}_snps_autosome_cleaned.bed",
        data_cache / "production_outputs/HWP/EUR_subjects.bed",
    )

    assert file_hashes_equal(
        tmp_path
        / f"population_level/EUR/controls_unrelated{pi}_maf{maf}_snps_autosome_cleaned.bim",
        data_cache / "production_outputs/HWP/EUR_subjects.bim",
    )

    assert file_hashes_equal(
        tmp_path
        / f"population_level/EUR/controls_unrelated{pi}_maf{maf}_snps_autosome_cleaned.fam",
        data_cache / "production_outputs/HWP/EUR_subjects.fam",
    )


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
@pytest.mark.slow
def test_phony_population_controls(tmp_path, conda_envs, qc_summary):
    # GIVEN: real data config, qc_summary table, all of the inputs to generate
    # the summary table, the European control list, and subject level plink data sets.
    # NOTE: I had to include the qc_summary inputs to get snakemake to run the
    # check point correctly.
    # NOTE: I am using the legacy EUR control list b/c of different sort orders.
    conda_envs.copy_env("plink2", tmp_path)
    conda_envs.copy_env("eigensoft", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy_sample_sheet()
        .add_user_files(entry_point="gtc", copy=False)
        .copy(
            "production_outputs/split_by_pop/EUR_subjects.bed", "population_level/EUR/subjects.bed",
        )
        .copy(
            "production_outputs/split_by_pop/EUR_subjects.bim", "population_level/EUR/subjects.bim",
        )
        .copy(
            "production_outputs/split_by_pop/EUR_subjects.fam", "population_level/EUR/subjects.fam",
        )
        .make_config(workflow_params={"subject_id_to_use": "PI_Subject_ID"})
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("common.smk")
            include: cfg.modules("plink_stats.smk")
            include: cfg.modules("plink_filters.smk")
            include: cfg.modules("subject_level_summary.smk")
            include: cfg.modules("population_level_summary.smk")

            rule all:
                input:
                    "population_level/controls.done",

            """
        )
    )
    (tmp_path / "sample_level").mkdir()
    shutil.copyfile(qc_summary, tmp_path / "sample_level/qc_summary.csv")

    with chdir(tmp_path):
        cfg = load_config()
        maf = cfg.config.software_params.maf_for_hwe
        pi = cfg.config.software_params.pi_hat_threshold

    # WHEN: run snakemake to get all population level and all population-control level results
    run_snakemake(tmp_path)

    # THEN:
    assert file_hashes_equal(
        tmp_path
        / f"population_level/EUR/controls_unrelated{pi}_maf{maf}_snps_autosome_cleaned.bed",
        data_cache / "production_outputs/HWP/EUR_subjects.bed",
    )
    assert file_hashes_equal(
        tmp_path
        / f"population_level/EUR/controls_unrelated{pi}_maf{maf}_snps_autosome_cleaned.hwe",
        data_cache / "production_outputs/HWP/EUR_subjects_qc.hwe",
    )
