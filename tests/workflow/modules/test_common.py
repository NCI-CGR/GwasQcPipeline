from textwrap import dedent

import pytest

from cgr_gwas_qc.testing import file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import FakeData, RealData


def test_bed_to_ped(tmp_path, conda_envs):
    # GIVEN:
    conda_envs.copy_env("plink2", tmp_path)
    (
        FakeData(tmp_path)
        .make_cgr_sample_sheet()
        .copy("plink/samples.bed", tmp_path / "samples.bed")
        .copy("plink/samples.bim", tmp_path / "samples.bim")
        .copy("plink/samples.fam", tmp_path / "samples.fam")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("common.smk")

            rule all:
                input:
                    "samples.ped",
                    "samples.map",
            """
        )
    )
    # WHEN:
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: the files exists
    assert (tmp_path / "samples.ped").exists()
    assert (tmp_path / "samples.map").exists()

    # NOTE: I tried comparing to the original PED/MAPs but plink is changing
    # marker order and allele order. This would require building a obs and exp
    # matrix using both PED/MAP. I assume PLINK works and testing the files
    # exist is good enough for now.


def test_eigenstrat_config(tmp_path):
    # GIVEN: a fake PED/MAP and config
    (
        FakeData(tmp_path)
        .make_cgr_sample_sheet()
        .copy("plink/samples.ped", tmp_path / "samples.ped")
        .copy("plink/samples.map", tmp_path / "samples.map")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("common.smk")

            rule all:
                input:
                    "samples.convert.par",
            """
        )
    )

    # WHEN: I run snakemake to create the par file
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: The config file is generated as expected
    obs_ = (tmp_path / "samples.convert.par").read_text()
    exp_ = dedent(
        """\
        genotypename: samples.ped
        snpname: samples.map
        indivname: samples.ped
        outputformat: EIGENSTRAT
        genooutfilename: samples.gen
        snpoutfilename: samples.snp
        indoutfilename: samples.ind
        familynames: NO
        """
    )

    assert obs_ == exp_


@pytest.mark.real_data
def test_eigensoft_convert(tmp_path, conda_envs):
    # GIVEN:
    conda_envs.copy_env("eigensoft", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .make_cgr_sample_sheet()
        .copy("production_outputs/pca/EUR_subjects_ld_prune.ped", tmp_path / "subjects.ped")
        .copy("production_outputs/pca/EUR_subjects_ld_prune.map", tmp_path / "subjects.map")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("common.smk")

            rule all:
                input:
                    "subjects.gen",
                    "subjects.snp",
                    "subjects.ind",
            """
        )
    )
    # WHEN:
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: the files exists
    file_hashes_equal(
        tmp_path / "subjects.gen",
        data_cache / "production_outputs/pca/EUR_subjects.eigenstratgeno",
    )

    file_hashes_equal(
        tmp_path / "subjects.snp", data_cache / "production_outputs/pca/EUR_subjects.snp"
    )

    file_hashes_equal(
        tmp_path / "subjects.ind", data_cache / "production_outputs/pca/EUR_subjects.ind"
    )


@pytest.mark.real_data
def test_eigensoft_smartpca(tmp_path, conda_envs):
    # GIVEN:
    conda_envs.copy_env("eigensoft", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .make_cgr_sample_sheet()
        .copy("production_outputs/pca/EUR_subjects.eigenstratgeno", tmp_path / "subjects.gen")
        .copy("production_outputs/pca/EUR_subjects.snp", tmp_path / "subjects.snp")
        .copy("production_outputs/pca/EUR_subjects.ind", tmp_path / "subjects.ind")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("common.smk")

            rule all:
                input:
                    "subjects.eigenvec",
            """
        )
    )
    # WHEN:
    run_snakemake(tmp_path)

    # THEN: the files exists
    file_hashes_equal(
        tmp_path / "subjects.eigenvec", data_cache / "production_outputs/pca/EUR_subjects.eigenvec",
    )
