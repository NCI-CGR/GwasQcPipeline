import pytest

from cgr_gwas_qc.testing import file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.real_data
@pytest.fixture(scope="session")
def call_rate_2_stats(tmp_path_factory, conda_envs):
    # GIVEN: A real data cache with samples after call rate 2 filtering.
    tmp_path = tmp_path_factory.mktemp("call_rate_2_stats")
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .add_sample_sheet()
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
