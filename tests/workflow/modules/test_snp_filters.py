import re
from typing import Tuple

import pytest

from cgr_gwas_qc.testing import file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import RealData


################################################################################
# MAF Filtering
################################################################################
def parse_plink_log_for_maf_counts(file_name) -> int:
    m = re.findall(r"\n(\d+) variants removed due to minor allele threshold", file_name.read_text())
    return int(m[0])


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_maf_filter(tmp_path, conda_envs):
    # GIVEN: A real data repository and a plink2 conda environmet
    conda_envs.copy_env("plink2", tmp_path)
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet()
        .copy("production_outputs/plink_filter_call_rate_2", "plink_filter_call_rate_2")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("snp_filters.smk")

            rule all:
                input: expand("snp_filters/maf/samples.{ext}", ext=["bed", "bim", "fam"])
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path)

    # THEN: The number of MAF filtered snps should be the same as the legacy run
    # Note: the legacy workflow does multiple filtering steps at the same time,
    # so I need to parse the log and just pull out the MAF filtering.
    obs_maf_count = parse_plink_log_for_maf_counts(tmp_path / "snp_filters/maf/samples.log")
    exp_maf_count = parse_plink_log_for_maf_counts(
        data_store / "production_outputs/ld_prune/ldPruneList.log"
    )
    assert obs_maf_count == exp_maf_count


################################################################################
# LD Pruning
################################################################################
def parse_plink_log_for_ld_pruning(file_name) -> Tuple[int, int]:
    m = re.findall(r"Pruning complete.\s+(\d+) of (\d+) variants removed.", file_name.read_text())
    pruned, kept = m[0]
    return int(pruned), int(kept)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_approx_ld_estimate(tmp_path, conda_envs):
    # GIVEN: A real data repository and a plink2 conda environmet
    # and running ld pruning
    conda_envs.copy_env("plink2", tmp_path)
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet()
        .copy("production_outputs/plink_filter_call_rate_2", "plink_filter_call_rate_2")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("snp_filters.smk")

            rule all:
                input: expand("snp_filters/ld_prune/ldPruneList.{ext}", ext=["prune.in", "prune.out"])
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path)

    # THEN: I expect the number of pruned snps to be the same.
    obs_pruned, _ = parse_plink_log_for_ld_pruning(
        tmp_path / "snp_filters/ld_prune/ldPruneList.log"
    )
    exp_pruned, _ = parse_plink_log_for_ld_pruning(
        data_store / "production_outputs/ld_prune/ldPruneList.log"
    )
    assert obs_pruned == exp_pruned


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_extract_ld_prune(tmp_path, conda_envs):
    # GIVEN:
    conda_envs.copy_env("plink2", tmp_path)
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet()
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("snp_filters.smk")

            rule all:
                input: expand("snp_filters/ld_prune/samples.{ext}", ext=["bed", "bim", "fam"])
            """
        )
    )
    (tmp_path / "snp_filters").mkdir(exist_ok=True)
    data_store.copy("production_outputs/plink_filter_call_rate_2", "snp_filters/maf")

    (tmp_path / "snp_filters/ld_prune").mkdir(exist_ok=True)
    data_store.copy(
        "production_outputs/ld_prune/ldPruneList.prune.in",
        "snp_filters/ld_prune/ldPruneList.prune.in",
    )
    data_store.copy(
        "production_outputs/ld_prune/ldPruneList.prune.out",
        "snp_filters/ld_prune/ldPruneList.prune.out",
    )
    data_store.copy(
        "production_outputs/ld_prune/ldPruneList.nosex", "snp_filters/ld_prune/ldPruneList.nosex"
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path)

    # THEN:
    obs_bed = tmp_path / "snp_filters/ld_prune/samples.bed"
    exp_bed = data_store / "production_outputs/ld_prune/samples.bed"
    assert file_hashes_equal(obs_bed, exp_bed)
