import re

import pytest

from cgr_gwas_qc import load_config
from cgr_gwas_qc.testing import chdir, file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import RealData


################################################################################
# Call Rate Filters
################################################################################
@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_call_rate_filters(tmp_path, conda_envs):
    """Check call rate filters."""
    # GIVEN: A real data repo and plink2 conda env
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy("production_outputs/plink_start/samples.bed", "sample_level/samples.bed")
        .copy("production_outputs/plink_start/samples.bim", "sample_level/samples.bim")
        .copy("production_outputs/plink_start/samples.fam", "sample_level/samples.fam")
        .make_config()
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("plink_filters.smk")

            rule all:
                input: expand("sample_level/call_rate_2/samples.{ext}", ext=["bed", "bim", "fam"])
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: The final file after filtering should match the file from the legacy workflow.
    obs_cr1 = tmp_path / "sample_level/call_rate_1/samples.bed"
    exp_cr1 = data_cache / "production_outputs/plink_filter_call_rate_1/samples.bed"
    assert file_hashes_equal(obs_cr1, exp_cr1)

    obs_cr2 = tmp_path / "sample_level/call_rate_2/samples.bed"
    exp_cr2 = data_cache / "production_outputs/plink_filter_call_rate_2/samples.bed"
    assert file_hashes_equal(obs_cr2, exp_cr2)


################################################################################
# MAF and LD Filters
################################################################################
@pytest.fixture(scope="session")
def maf_and_ld_outputs(tmp_path_factory, conda_envs):
    tmp_path = tmp_path_factory.mktemp("maf_and_ld_outputs")

    # GIVEN: A real data repository and a plink2 conda environmet
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples.bed",
            "sample_level/call_rate_2/samples.bed",
        )
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples.bim",
            "sample_level/call_rate_2/samples.bim",
        )
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples.fam",
            "sample_level/call_rate_2/samples.fam",
        )
        .make_config()
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("plink_filters.smk")

            rule all:
                input:
                    "sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.bed".format(
                        maf=cfg.config.software_params.maf_for_ibd,
                        ld=cfg.config.software_params.ld_prune_r2
                    )
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path, keep_temp=True)

    return data_cache, tmp_path


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_maf_filter(maf_and_ld_outputs):
    # GIVEN: Outputs from the maf_filter rule
    data_cache, tmp_path = maf_and_ld_outputs

    # THEN: The number of MAF filtered snps should be the same as the legacy run
    # Note: the legacy workflow does multiple filtering steps at the same time,
    # so I need to parse the log and just pull out the MAF filtering.
    def parse_plink_log_for_maf_counts(file_name) -> int:
        m = re.findall(
            r"\n(\d+) variants removed due to minor allele threshold", file_name.read_text()
        )
        return int(m[0])

    with chdir(tmp_path):
        cfg = load_config(pytest=True)
        maf = cfg.config.software_params.maf_for_ibd

    obs_maf_count = parse_plink_log_for_maf_counts(
        tmp_path / f"sample_level/call_rate_2/samples_maf{maf}.log"
    )
    exp_maf_count = parse_plink_log_for_maf_counts(
        data_cache / "production_outputs/ld_prune/ldPruneList.log"
    )
    assert obs_maf_count == exp_maf_count


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_approx_ld(maf_and_ld_outputs):
    # GIVEN: A real data repository and a plink2 conda environmet
    data_cache, tmp_path = maf_and_ld_outputs

    # THEN: I expect the number of pruned snps to be the same.
    def parse_plink_log_for_ld_pruning(file_name):
        m = re.findall(
            r"Pruning complete.\s+(\d+) of (\d+) variants removed.", file_name.read_text()
        )
        pruned, kept = m[0]
        return int(pruned), int(kept)

    with chdir(tmp_path):
        cfg = load_config(pytest=True)
        maf, ld = cfg.config.software_params.maf_for_ibd, cfg.config.software_params.ld_prune_r2

    obs_pruned, _ = parse_plink_log_for_ld_pruning(
        tmp_path / f"sample_level/call_rate_2/samples_maf{maf}_ld{ld}.log"
    )
    exp_pruned, _ = parse_plink_log_for_ld_pruning(
        data_cache / "production_outputs/ld_prune/ldPruneList.log"
    )
    assert obs_pruned == exp_pruned


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_ld_prune(maf_and_ld_outputs):
    # GIVEN:
    data_cache, tmp_path = maf_and_ld_outputs

    # THEN:
    with chdir(tmp_path):
        cfg = load_config(pytest=True)
        maf, ld = cfg.config.software_params.maf_for_ibd, cfg.config.software_params.ld_prune_r2

    obs_bed = tmp_path / f"sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.bed"
    exp_bed = data_cache / "production_outputs/ld_prune/samples.bed"
    assert file_hashes_equal(obs_bed, exp_bed)
