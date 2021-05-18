import shutil

import pytest

from cgr_gwas_qc.testing import file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.workflow
@pytest.mark.regression
@pytest.mark.real_data
def test_create_subjects(tmp_path, conda_envs, sample_qc_csv, sample_concordance_csv):
    # GIVEN:
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

            include: cfg.modules("subject_level_summary")

            rule all:
                input:
                    "subject_level/subjects.bed",
                    "subject_level/subjects.bim",
                    "subject_level/subjects.fam",
            """
        )
    )

    shutil.copyfile(sample_qc_csv, tmp_path / "sample_level/sample_qc.csv")
    (tmp_path / "sample_level/concordance").mkdir()
    shutil.copyfile(sample_concordance_csv, tmp_path / "sample_level/concordance/summary.csv")

    # WHEN:
    run_snakemake(tmp_path)

    # THEN:
    file_hashes_equal(
        tmp_path / "subject_level/subjects.bed",
        data_cache / "production_outputs/subject_level/subjects.bed",
    )
    file_hashes_equal(
        tmp_path / "subject_level/subjects.bim",
        data_cache / "production_outputs/subject_level/subjects.bim",
    )
    file_hashes_equal(
        tmp_path / "subject_level/subjects.fam",
        data_cache / "production_outputs/subject_level/subjects.fam",
    )


@pytest.mark.workflow
@pytest.mark.real_data
def test_remove_related_subjects(tmp_path, conda_envs):
    # GIVEN: Real subject level data
    conda_envs.copy_env("plink2", tmp_path)
    (
        RealData(tmp_path)
        .copy(
            "production_outputs/subject_level/subjects.bed",
            "subject_level/subjects.bed",
        )
        .copy(
            "production_outputs/subject_level/subjects.bim",
            "subject_level/subjects.bim",
        )
        .copy(
            "production_outputs/subject_level/subjects.fam",
            "subject_level/subjects.fam",
        )
        .make_config(software_params={"pi_hat_threshold": 0.16})  # ensure some subjects are related
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("plink_filters")
            include: cfg.modules("plink_stats")
            include: cfg.modules("subject_level_summary")

            rule all:
                input:
                    expand(
                        "subject_level/subjects_unrelated.{ext}",
                        pi=cfg.config.software_params.pi_hat_threshold,
                        ext=["bed", "bim", "fam"]
                    )
            """
        )
    )

    # WHEN: Prune related subjects.
    run_snakemake(tmp_path)

    # THEN: Given my test data and the `pi_hat_threshold = 0.16` I expect 7 subjects to be remove.
    assert (
        len(
            (tmp_path / "subject_level/related_subjects_to_remove.txt")
            .read_text()
            .strip()
            .splitlines()
        )
        == 6
    )

    # The plink log should say I 177 subjects remaining. (184 - 7 = 177)
    log = (tmp_path / "subject_level/subjects_unrelated.log").read_text()
    assert "178 people remaining." in log
