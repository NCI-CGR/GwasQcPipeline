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
def files_for_upload(tmp_path_factory, sample_qc_csv):
    tmp_path = tmp_path_factory.mktemp("files_for_upload")
    data_cache = (
        RealData(tmp_path)
        .copy("original_data/manifest_full.csv", "SR001-001_00_AnalysisManifest_0000000.csv")
        .copy(
            "production_outputs/concordance/KnownReplicates.csv",
            "sample_level/concordance/KnownReplicates.csv",
        )
        .copy(
            "production_outputs/concordance/UnknownReplicates.csv",
            "sample_level/concordance/UnknownReplicates.csv",
        )
        .make_cgr_sample_sheet()
        .make_config(
            sample_sheet="SR001-001_00_AnalysisManifest_0000000.csv",
            user_files=dict(output_pattern="{prefix}/SR001-001_00_{file_type}_0000000.{ext}"),
        )
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
    shutil.copyfile(sample_qc_csv, tmp_path / "sample_level/sample_qc.csv")

    # WHEN:
    run_snakemake(tmp_path)

    return data_cache, tmp_path


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
    obs_df = (
        pd.read_csv(tmp_path / "files_for_lab/SR001-001_00_LimsUpload_0000000.csv")
        .set_index("Sample_ID")
        .fillna({"Expected Replicate Discordance": False, "Unexpected Replicate": False})
    )
    exp_df = (
        pd.read_csv(
            data_cache
            / "production_outputs/files_for_lab/SR0446-001_12_LimsUpload_1011201995419_casecontrol_20191011.csv"
        )
        .rename({"Sample ID": "Sample_ID", "Call Rate": "Call_Rate_Initial"}, axis=1)
        .set_index("Sample_ID")
    )

    assert_frame_equal(exp_df, obs_df, check_dtype=False, check_like=True)


@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.real_data
def test_identifiler_needed(files_for_upload):
    # GIVEN: Real data after running the delivery module
    data_cache, tmp_path = files_for_upload

    # THEN: The Identifiler table for upload should be identical with production outputs.
    obs_ = pd.read_csv(tmp_path / "files_for_lab/SR001-001_00_Identifiler_0000000.csv")
    exp_ = pd.read_csv(
        data_cache
        / "production_outputs/files_for_lab/SR0446-001_12_Identifiler_1011201995419_casecontrol_20191011.csv"
    ).replace({"Sex Discordant": "Sex Discordance"})

    assert_frame_equal(obs_, exp_, check_dtype=False)


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
def files_for_deliver(tmp_path_factory, sample_qc_csv):
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
        .make_cgr_sample_sheet()
        .make_config(
            sample_sheet="SR001-001_00_AnalysisManifest_0000000.csv",
            user_files=dict(output_pattern="{prefix}/SR001-001_00_{file_type}_0000000.{ext}"),
        )
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
    shutil.copyfile(sample_qc_csv, tmp_path / "sample_level/sample_qc.csv")

    # HWE results for Europeans b/c they are the only non-empty files
    with chdir(tmp_path):
        cfg = load_config(pytest=True)
        maf = cfg.config.software_params.maf_for_hwe

    data_cache.copy(
        "production_outputs/HWP/EUR_subjects_qc.hwe",
        f"population_level/EUR/controls_maf{maf}_snps_autosome_cleaned.hwe",
    )
    (tmp_path / "population_level/per_population_controls_qc.done").write_text(
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
        cfg = load_config(pytest=True)
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


################################################################################
# Reports
################################################################################
@pytest.mark.real_data
def test_qc_report_table(
    sample_concordance_csv,
    sample_qc_csv,
    subject_qc_csv,
    agg_population_concordance_csv,
    population_qc_csv,
    tmp_path,
):
    (
        RealData(tmp_path)
        .make_cgr_sample_sheet()
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("delivery.smk")

            rule all:
                input:
                    "deliver/test_QC_Report_test.xlsx"
            """
        )
    )
    (tmp_path / "sample_level/ancestry").mkdir(parents=True)
    (tmp_path / "sample_level/concordance").mkdir(parents=True)
    (tmp_path / "subject_level").mkdir(parents=True)
    (tmp_path / "population_level").mkdir(parents=True)

    shutil.copy(sample_concordance_csv, tmp_path / "sample_level/concordance/summary.csv")
    (tmp_path / "sample_level/ancestry/graf_populations.txt").write_text(
        "DS No.\tSample\t#SNPs\tGD1\tGD2\tGD3\tGD4\tF(%)\tE(%)\tA(%)\tAfrican\tEuropean\tAsian\tMexican\tIndian-Pakistani\n"
    )
    shutil.copy(sample_qc_csv, tmp_path / "sample_level/sample_qc.csv")
    shutil.copy(subject_qc_csv, tmp_path / "subject_level/subject_qc.csv")
    shutil.copy(population_qc_csv, tmp_path / "population_level/population_qc.csv")
    shutil.copy(agg_population_concordance_csv, tmp_path / "population_level/concordance.csv")

    run_snakemake(tmp_path)
