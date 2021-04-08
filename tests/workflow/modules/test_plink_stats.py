import pytest

from cgr_gwas_qc import load_config
from cgr_gwas_qc.testing import chdir, file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.real_data
@pytest.mark.workflow
@pytest.fixture(scope="session")
def call_rate_2_stats(tmp_path_factory, conda_envs):
    # GIVEN: A real data cache with samples after call rate 2 filtering.
    tmp_path = tmp_path_factory.mktemp("call_rate_2_stats")
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

            include: cfg.modules("plink_stats.smk")

            rule all:
                input:
                    "sample_level/call_rate_2/samples.imiss",
                    "sample_level/call_rate_2/samples.lmiss",
                    "sample_level/call_rate_2/samples.sexcheck",
                    "sample_level/call_rate_2/samples.frq",
                    "sample_level/call_rate_2/samples.hwe",
            """
        )
    )
    # WHEN: I run snakemake to create all stats
    run_snakemake(tmp_path)

    return data_cache, tmp_path


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_plink_stats_call_rate(call_rate_2_stats):
    # GIVEN: outputs from running all plink stats on call rate 2 samples
    data_cache, tmp_path = call_rate_2_stats

    # THEN: The outputs should match the legacy workflow
    assert file_hashes_equal(
        tmp_path / "sample_level/call_rate_2/samples.imiss",
        data_cache / "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss",
    )
    assert file_hashes_equal(
        tmp_path / "sample_level/call_rate_2/samples.lmiss",
        data_cache / "production_outputs/plink_filter_call_rate_2/samples_filter2.lmiss",
    )


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_plink_stats_sexcheck(call_rate_2_stats):
    # GIVEN: outputs from running all plink stats on call rate 2 samples
    data_cache, tmp_path = call_rate_2_stats

    # THEN: The outputs should match the legacy workflow
    assert file_hashes_equal(
        tmp_path / "sample_level/call_rate_2/samples.sexcheck",
        data_cache / "production_outputs/plink_filter_call_rate_2/samples_filter2.sexcheck",
    )


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_plink_stats_allele_freq(call_rate_2_stats):
    # GIVEN: outputs from running all plink stats on call rate 2 samples
    data_cache, tmp_path = call_rate_2_stats

    # THEN: The outputs should match the legacy workflow
    assert file_hashes_equal(
        tmp_path / "sample_level/call_rate_2/samples.frq",
        data_cache / "production_outputs/plink_filter_call_rate_2/samples_filter2.frq",
    )


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_plink_stats_hardy(call_rate_2_stats):
    # GIVEN: outputs from running all plink stats on call rate 2 samples
    data_cache, tmp_path = call_rate_2_stats

    # THEN: The outputs should match the legacy workflow
    assert file_hashes_equal(
        tmp_path / "sample_level/call_rate_2/samples.hwe",
        data_cache / "production_outputs/plink_filter_call_rate_2/samples_filter2.hwe",
    )


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_plink_stats_ibd(tmp_path, conda_envs):
    """IBD estimation

    Unlike the other stats, IBD estimation needs to be done on LD-pruned
    outputs.
    """
    # GIVEN: plink2, and real data and config
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .make_config()
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("plink_stats.smk")
            include: cfg.modules("plink_filters.smk")

            rule all:
                input:
                    "sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.genome".format(
                        maf=cfg.config.software_params.maf_for_ibd,
                        ld=cfg.config.software_params.ld_prune_r2
                    )
            """
        )
    )
    with chdir(tmp_path):
        cfg = load_config(pytest=True)
        maf, ld = cfg.config.software_params.maf_for_ibd, cfg.config.software_params.ld_prune_r2

    # And outputs from LD pruning
    (
        data_cache.copy(
            "production_outputs/ld_prune/samples.bed",
            f"sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.bed",
        )
        .copy(
            "production_outputs/ld_prune/samples.bim",
            f"sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.bim",
        )
        .copy(
            "production_outputs/ld_prune/samples.fam",
            f"sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.fam",
        )
    )

    # WHEN: I run snakemake to generate IBD outputs
    run_snakemake(tmp_path)

    # THEN: The legacy and new outputs should match exactly
    obs_ = tmp_path / f"sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.genome"
    exp_ = data_cache / "production_outputs/ibd/samples.genome"
    assert file_hashes_equal(obs_, exp_)
