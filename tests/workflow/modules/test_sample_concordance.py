import pytest

from cgr_gwas_qc.testing import file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_plink_ibd(tmp_path, conda_envs):
    # GIVEN:
    conda_envs.copy_env("plink2", tmp_path)
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .make_config()
        .copy("production_outputs/ld_prune", "snp_filters/ld_prune")
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("sample_concordance.smk")

            rule all:
                input:
                    expand("sample_concordance/ibd/samples.{ext}", ext=["genome", "nosex"])
            """
        )
    )

    # WHEN:
    run_snakemake(tmp_path)

    # THEN:
    obs_file = tmp_path / "sample_concordance/ibd/samples.genome"
    exp_file = data_store / "production_outputs/ibd/samples.genome"
    assert file_hashes_equal(obs_file, exp_file)


@pytest.mark.real_data
@pytest.mark.workflow
def test_sample_concordance(tmp_path):
    # GIVEN: Real data
    (
        RealData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .make_config()
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss",
            "plink_filter_call_rate_2/samples.imiss",
        )
        .copy("production_outputs/ibd/samples.genome", "sample_concordance/ibd/samples.genome")
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("sample_concordance.smk")

            rule all:
                input:
                    "sample_concordance/KnownReplicates.csv",
                    "sample_concordance/InternalQcKnown.csv",
                    "sample_concordance/StudySampleKnown.csv",
                    "sample_concordance/UnknownReplicates.csv",
            """
        )
    )

    # WHEN: we run snakemake looking for the replicate output files.
    run_snakemake(tmp_path)

    # THEN: these files are created.
    # Note: regression tests are executed as part of the
    # `known_concorant_samples.py` testing
    assert (tmp_path / "sample_concordance/KnownReplicates.csv").exists()
    assert (tmp_path / "sample_concordance/InternalQcKnown.csv").exists()
    assert (tmp_path / "sample_concordance/StudySampleKnown.csv").exists()
    assert (tmp_path / "sample_concordance/UnknownReplicates.csv").exists()
