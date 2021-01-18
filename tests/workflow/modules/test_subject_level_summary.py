import pytest

from cgr_gwas_qc.testing import file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.workflow
@pytest.mark.regression
@pytest.mark.real_data
def test_create_subjects(tmp_path, conda_envs):
    # GIVEN:
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .add_sample_sheet()
        .copy(
            "production_outputs/subject_level/SampleUsedforSubject.csv",
            "sample_level/samples_used_for_subjects.csv",
        )
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

            include: cfg.modules("subject_level_summary.smk")

            rule all:
                input:
                    "subject_level/subjects.bed",
                    "subject_level/subjects.bim",
                    "subject_level/subjects.fam",
            """
        )
    )

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
