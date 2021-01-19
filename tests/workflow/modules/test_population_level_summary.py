import shutil

import pytest

from cgr_gwas_qc.testing import file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_subjects_per_population(tmp_path, qc_summary):
    # GIVEN
    data_cache = (
        RealData(tmp_path)
        .add_sample_sheet()
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("population_level_summary.smk")

            rule all:
                input:
                    "population_level/AFR/subject_list.txt",
                    "population_level/ASN/subject_list.txt",
                    "population_level/EUR/subject_list.txt",
            """
        )
    )
    (tmp_path / "sample_level").mkdir()
    shutil.copyfile(qc_summary, tmp_path / "sample_level/qc_summary.csv")

    # WHEN
    run_snakemake(tmp_path)

    # THEN
    file_hashes_equal(
        tmp_path / "population_level/AFR/subject_list.txt",
        data_cache / "production_outputs/split_by_pop/AFR.keep.txt",
    )

    file_hashes_equal(
        tmp_path / "population_level/ASN/subject_list.txt",
        data_cache / "production_outputs/split_by_pop/ASN.keep.txt",
    )

    file_hashes_equal(
        tmp_path / "population_level/EUR/subject_list.txt",
        data_cache / "production_outputs/split_by_pop/EUR.keep.txt",
    )


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_plink_split_population(tmp_path, conda_envs, qc_summary):
    # GIVEN
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .add_sample_sheet()
        .copy(
            "production_outputs/split_by_pop/AFR.keep.txt", "population_level/AFR/subject_list.txt",
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
                    "population_level/AFR/subjects.bed",
                    "population_level/AFR/subjects.bim",
                    "population_level/AFR/subjects.fam",
            """
        )
    )
    (tmp_path / "sample_level").mkdir()
    shutil.copyfile(qc_summary, tmp_path / "sample_level/qc_summary.csv")

    # WHEN
    run_snakemake(tmp_path)

    # THEN
    file_hashes_equal(
        tmp_path / "population_level/AFR/subjects.bed",
        data_cache / "production_outputs/split_by_pop/AFR_subjects.bed",
    )

    file_hashes_equal(
        tmp_path / "population_level/AFR/subjects.bim",
        data_cache / "production_outputs/split_by_pop/AFR_subjects.bim",
    )

    file_hashes_equal(
        tmp_path / "population_level/AFR/subjects.fam",
        data_cache / "production_outputs/split_by_pop/AFR_subjects.fam",
    )
