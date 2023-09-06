import pytest

from cgr_gwas_qc.testing import run_snakemake
from cgr_gwas_qc.testing.comparison import file_hashes_equal
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.workflow
@pytest.mark.regression
@pytest.mark.real_data
def test_eigensoft_convert(tmp_path, conda_envs):
    # GIVEN:
    conda_envs.copy_env("eigensoft", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .make_cgr_sample_sheet()
        .make_config()
        .copy("legacy_outputs/pca/EUR_subjects_ld_prune.ped", "samples.ped")
        .copy("legacy_outputs/pca/EUR_subjects_ld_prune.map", "samples.map")
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            rule all:
                input:
                    "samples.convert.par",
                    "samples.gen",
                    "samples.snp",
                    "samples.ind",

            module eigensoft:
                snakefile: cfg.modules("eigensoft")

            use rule convert from eigensoft

            """
        )
    )
    # WHEN:
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: the files exists
    file_hashes_equal(
        data_cache / "legacy_outputs/pca/EUR_subjects.convertEigen.par",
        tmp_path / "samples.convert.par",
    )
    file_hashes_equal(
        data_cache / "legacy_outputs/pca/EUR_subjects.eigenstratgeno", tmp_path / "samples.gen"
    )
    file_hashes_equal(data_cache / "legacy_outputs/pca/EUR_subjects.snp", tmp_path / "samples.snp")
    file_hashes_equal(data_cache / "legacy_outputs/pca/EUR_subjects.ind", tmp_path / "samples.ind")


@pytest.mark.workflow
@pytest.mark.regression
@pytest.mark.real_data
def test_eigensoft_smartpca(tmp_path, conda_envs):
    # GIVEN:
    conda_envs.copy_env("eigensoft", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .make_cgr_sample_sheet()
        .copy("legacy_outputs/pca/EUR_subjects.eigenstratgeno", tmp_path / "samples.gen")
        .copy("legacy_outputs/pca/EUR_subjects.snp", tmp_path / "samples.snp")
        .copy("legacy_outputs/pca/EUR_subjects.ind", tmp_path / "samples.ind")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            rule all:
                input:
                    "samples.pca.par",
                    "samples.eigenvec",

            module eigensoft:
                snakefile: cfg.modules("eigensoft")

            use rule smartpca from eigensoft

            """
        )
    )
    # WHEN:
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: the files exists
    file_hashes_equal(
        data_cache / "legacy_outputs/pca/EUR_subjects.smartpca.par", tmp_path / "samples.pca.par",
    )
    file_hashes_equal(
        data_cache / "legacy_outputs/pca/EUR_subjects.eigenvec", tmp_path / "samples.eigenvec",
    )
