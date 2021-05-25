import pytest

from cgr_gwas_qc.testing import run_snakemake
from cgr_gwas_qc.testing.comparison import file_hashes_equal
from cgr_gwas_qc.testing.data import RealData

################################################################################
# Legacy Regression Tests
################################################################################


################################################################################
# Workflow Tests
################################################################################
@pytest.mark.workflow
@pytest.mark.real_data
@pytest.fixture(scope="module")
def delivery(pytestconfig, tmp_path_factory, conda_envs):
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    tmp_path = tmp_path_factory.mktemp("delivery")
    conda_envs.copy_env("pandoc", tmp_path)
    (
        RealData(tmp_path, full_sample_sheet=False)
        .copy("dev_outputs/cgr_sample_sheet.csv", "cgr_sample_sheet.csv")
        .copy("dev_outputs/config.yml", "config.yml")
        .copy("dev_outputs/sample_level", "sample_level")
        .copy("dev_outputs/subject_level", "subject_level")
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            module delivery:
                snakefile: cfg.subworkflow("delivery")
                config: {}

            use rule * from delivery
            """
        )
    )

    run_snakemake(tmp_path, keep_temp=True)

    return tmp_path


@pytest.mark.parametrize(
    "filename",
    [
        "files_for_lab/all_sample_qc.csv",
        "files_for_lab/Identifiler.csv",
        "files_for_lab/KnownReplicates.csv",
        "files_for_lab/LimsUpload.csv",
        "files_for_lab/UnknownReplicates.csv",
        "delivery/samples.bed",
        "delivery/samples.bim",
        "delivery/samples.fam",
        "delivery/SampleUsedforSubject.csv",
        "delivery/subjects.bed",
        "delivery/subjects.bim",
        "delivery/subjects.fam",
        pytest.param(
            "delivery/qc_report.md",
            marks=pytest.mark.xfail(reason="contains dates and paths that change"),
        ),
        pytest.param(
            "delivery/QC_Report.docx",
            marks=pytest.mark.xfail(reason="contains dates and paths that change"),
        ),
        pytest.param(
            "delivery/QC_Report.xlsx",
            marks=pytest.mark.xfail(reason="contains dates and paths that change"),
        ),
        pytest.param(
            "delivery/HWP.zip",
            marks=pytest.mark.xfail(reason="zip files contain header data that will change"),
        ),
    ],
)
@pytest.mark.workflow
@pytest.mark.real_data
def test_file_hashes_equal(real_data_cache, delivery, filename):
    assert file_hashes_equal(real_data_cache / "dev_outputs" / filename, delivery / filename)
