import shutil

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from cgr_gwas_qc import load_config
from cgr_gwas_qc.testing import chdir, file_hashes_equal, run_snakemake, sorted_file_equal
from cgr_gwas_qc.testing.data import RealData


################################################################################
# Files For Lab
################################################################################
@pytest.mark.workflow
@pytest.mark.real_data
@pytest.fixture(scope="session")
def files_for_upload(tmp_path_factory, qc_summary):
    tmp_path = tmp_path_factory.mktemp("files_for_upload")
    data_cache = (
        RealData(tmp_path)
        .copy("original_data/manifest_full.csv", "SR001-001_00_AnalysisManifest_0000000.csv",)
        .copy(
            "production_outputs/concordance/KnownReplicates.csv",
            "sample_level/concordance/KnownReplicates.csv",
        )
        .copy(
            "production_outputs/concordance/UnknownReplicates.csv",
            "sample_level/concordance/UnknownReplicates.csv",
        )
        .make_config(sample_sheet="SR001-001_00_AnalysisManifest_0000000.csv")
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("delivery.smk")

            rule all:
                input:
                    "files_for_lab/SR001-001_00_all_sample_qc_0000000.csv",
                    "files_for_lab/SR001-001_00_LimsUpload_0000000.csv",
                    "files_for_lab/SR001-001_00_Identifiler_0000000.csv",
                    "files_for_lab/SR001-001_00_KnownReplicates_0000000.csv",
                    "files_for_lab/SR001-001_00_UnknownReplicates_0000000.csv",
            """
        )
    )
    shutil.copyfile(qc_summary, tmp_path / "sample_level/qc_summary.csv")

    # WHEN:
    run_snakemake(tmp_path)

    return data_cache, tmp_path


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_lab_sample_level_qc_report(files_for_upload):
    # GIVEN: Real data after running the delivery module
    data_cache, tmp_path = files_for_upload

    # THEN: The QC table for upload should be identical with production outputs.
    # NOTE: After removing new columns
    obs_ = pd.read_csv(tmp_path / "files_for_lab/SR001-001_00_all_sample_qc_0000000.csv")
    exp_ = pd.read_csv(
        data_cache
        / "production_outputs/files_for_lab/SR0446-001_12_all_sample_qc_1011201995419_casecontrol_20191011.csv"
    )

    assert_frame_equal(
        exp_, obs_.reindex(exp_.columns, axis=1)  # Exclude columns not in the original table.
    )


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_lab_lims_upload(files_for_upload):
    # GIVEN: Real data after running the delivery module
    data_cache, tmp_path = files_for_upload

    # THEN: The LIMS table for upload should be identical with production outputs.
    # NOTE: The legacy workflow uses `Sample ID` instead of `Sample_ID`. The
    # dev workflow uses `Sample_ID` so I need to rename production outputs
    # prior to comparisons.
    obs_ = pd.read_csv(tmp_path / "files_for_lab/SR001-001_00_LimsUpload_0000000.csv")
    exp_ = pd.read_csv(
        data_cache
        / "production_outputs/files_for_lab/SR0446-001_12_LimsUpload_1011201995419_casecontrol_20191011.csv"
    ).rename({"Sample ID": "Sample_ID"}, axis=1)

    assert_frame_equal(obs_, exp_)


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_identifiler_needed(files_for_upload):
    # GIVEN: Real data after running the delivery module
    data_cache, tmp_path = files_for_upload

    # THEN: The Identifiler table for upload should be identical with production outputs.
    # NOTE: The legacy workflow uses `Identifiler Reason` where the dev workflow uses `Identifiler_Reason`.
    # I need to rename production outputs prior to comparisons.
    obs_ = pd.read_csv(tmp_path / "files_for_lab/SR001-001_00_Identifiler_0000000.csv")
    exp_ = pd.read_csv(
        data_cache
        / "production_outputs/files_for_lab/SR0446-001_12_Identifiler_1011201995419_casecontrol_20191011.csv"
    ).rename({"Identifiler Reason": "Identifiler_Reason"}, axis=1)

    assert_frame_equal(obs_, exp_)


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_known_replicate(files_for_upload):
    # GIVEN: Real data after running the delivery module
    data_cache, tmp_path = files_for_upload

    # THEN: The KnownReplicates table for upload should be identical with production outputs.
    assert file_hashes_equal(
        tmp_path / "files_for_lab/SR001-001_00_KnownReplicates_0000000.csv",
        data_cache
        / "production_outputs/files_for_lab/SR0446-001_12_KnownReplicates_1011201995419_casecontrol_20191011.csv",
    )


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_unknown_replicate(files_for_upload):
    # GIVEN: Real data after running the delivery module
    data_cache, tmp_path = files_for_upload

    # THEN: The KnownReplicates table for upload should be identical with production outputs.
    assert file_hashes_equal(
        tmp_path / "files_for_lab/SR001-001_00_UnknownReplicates_0000000.csv",
        data_cache
        / "production_outputs/files_for_lab/SR0446-001_12_UnknownReplicates_1011201995419_casecontrol_20191011.csv",
    )


################################################################################
# Deliverables
################################################################################
@pytest.mark.workflow
@pytest.mark.real_data
@pytest.fixture(scope="session")
def files_for_deliver(tmp_path_factory, qc_summary):
    tmp_path = tmp_path_factory.mktemp("files_for_deliver")

    # GIVEN: Real data
    data_cache = (
        RealData(tmp_path)
        # Sample files
        .copy("production_outputs/plink_start/samples.bed", "sample_level/samples.bed")
        .copy("production_outputs/plink_start/samples.bim", "sample_level/samples.bim")
        .copy("production_outputs/plink_start/samples.fam", "sample_level/samples.fam")
        # Subject files
        .copy("production_outputs/subject_level/samples.bed", "subject_level/samples.bed")
        .copy("production_outputs/subject_level/samples.bim", "subject_level/samples.bim")
        .copy("production_outputs/subject_level/samples.fam", "subject_level/samples.fam")
        # Use CGR manifest naming scheme
        .copy("original_data/manifest_full.csv", "SR001-001_00_AnalysisManifest_0000000.csv",)
        .make_config(sample_sheet="SR001-001_00_AnalysisManifest_0000000.csv")
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("delivery.smk")

            rule all:
                input:
                    "deliver/SR001-001_00_AnalysisManifest_0000000.csv",
                    "deliver/HWP.zip",
                    "deliver/samples.bed",
                    "deliver/samples.bim",
                    "deliver/samples.fam",
                    "deliver/subjects.bed",
                    "deliver/subjects.bim",
                    "deliver/subjects.fam",
                    "deliver/SampleUsedforSubject.csv",
            """
        )
    )
    # The QC Summary
    shutil.copyfile(qc_summary, tmp_path / "sample_level/qc_summary.csv")

    # HWE results for Europeans b/c they are the only non-empty files
    with chdir(tmp_path):
        cfg = load_config()
        maf = cfg.config.software_params.maf_for_hwe

    data_cache.copy(
        "production_outputs/HWP/EUR_subjects_qc.hwe",
        f"population_level/EUR/controls_maf{maf}_snps_autosome_cleaned.hwe",
    )
    (tmp_path / "population_level/controls.done").write_text(
        f"population_level/EUR/controls_maf{maf}_snps_autosome_cleaned.hwe"
    )

    # WHEN:
    run_snakemake(tmp_path)

    return data_cache, tmp_path


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_deliver_manifest(files_for_deliver):
    # GIVEN: Files in deliver folder
    data_cache, deliver_path = files_for_deliver

    # THEN: the manifest should match b/c it is just copied
    assert file_hashes_equal(
        data_cache / "original_data/manifest_full.csv",
        next(deliver_path.glob("*AnalysisManifest*")),
    )


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_deliver_hwp(files_for_deliver, tmp_path):
    from zipfile import ZipFile

    # GIVEN: Files in deliver folder
    data_cache, deliver_path = files_for_deliver

    # WHEN: Unzip `deliver/HWP.zip`
    with chdir(deliver_path):
        cfg = load_config()
        maf = cfg.config.software_params.maf_for_hwe

    with ZipFile(deliver_path / "deliver/HWP.zip") as myzip:
        myzip.extract(f"HWP/controls_maf{maf}_snps_autosome_cleaned.hwe", tmp_path)

    # THEN: The European HWE file should match b/c it is just copied
    assert file_hashes_equal(
        data_cache / "production_outputs/HWP/EUR_subjects_qc.hwe",
        tmp_path / f"HWP/controls_maf{maf}_snps_autosome_cleaned.hwe",
    )


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_deliver_original_sample_data(files_for_deliver):
    # GIVEN: Files in deliver folder
    data_cache, deliver_path = files_for_deliver

    # THEN: the sample data should match b/c they are just copied
    assert file_hashes_equal(
        data_cache / "production_outputs/plink_start/samples.bed",
        deliver_path / "deliver/samples.bed",
    )
    assert file_hashes_equal(
        data_cache / "production_outputs/plink_start/samples.bim",
        deliver_path / "deliver/samples.bim",
    )
    assert file_hashes_equal(
        data_cache / "production_outputs/plink_start/samples.fam",
        deliver_path / "deliver/samples.fam",
    )


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_deliver_subject_data(files_for_deliver):
    # GIVEN: Files in deliver folder
    data_cache, deliver_path = files_for_deliver

    # THEN: the subject data should match b/c they are just copied
    assert file_hashes_equal(
        data_cache / "production_outputs/subject_level/samples.bed",
        deliver_path / "deliver/subjects.bed",
    )
    assert file_hashes_equal(
        data_cache / "production_outputs/subject_level/samples.bim",
        deliver_path / "deliver/subjects.bim",
    )
    assert file_hashes_equal(
        data_cache / "production_outputs/subject_level/samples.fam",
        deliver_path / "deliver/subjects.fam",
    )


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_deliver_subject_list(files_for_deliver):
    # GIVEN: Files in deliver folder
    data_cache, deliver_path = files_for_deliver

    # THEN: the subject list should match but not necessarily in the same order.
    assert sorted_file_equal(
        data_cache / "production_outputs/subject_level/SampleUsedforSubject.csv",
        deliver_path / "deliver/SampleUsedforSubject.csv",
    )
