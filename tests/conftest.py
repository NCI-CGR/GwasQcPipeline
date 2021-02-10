from pathlib import Path

import pandas as pd
import pysam
import pytest

from cgr_gwas_qc.parsers.illumina import BeadPoolManifest, GenotypeCalls
from cgr_gwas_qc.parsers.sample_sheet import SampleSheet
from cgr_gwas_qc.testing.conda import CondaEnv
from cgr_gwas_qc.testing.data import FakeData, RealData


##################################################################################
# Pytest Configuration and Hooks
##################################################################################
def pytest_addoption(parser):
    """Add pytest command line options."""
    parser.addoption(
        "--real-data", action="store_true", default=False, help="Run tests with real data."
    )

    parser.addoption("--no-slow", action="store_true", default=False, help="Don't run slow tests.")

    parser.addoption(
        "--sync-real-data",
        action="store_true",
        default=False,
        help="Synchronize with real data store.",
    )


def pytest_collection_modifyitems(config, items):
    """Modify tests based on command line options.

    Based on an example from the pytest website:
        https://docs.pytest.org/en/6.0.1/example/simple.html#control-skipping-of-tests-according-to-command-line-option
    """
    if config.getoption("--sync-real-data"):
        RealData(sync=True)  # Make sure real data is downloaded

    if not config.getoption("--real-data"):
        # Skip test that are marked with `real_data` unless the command line
        # flag `--real-data` is provided
        skip_real_data = pytest.mark.skip(reason="need --real-data option to run")
        for item in items:
            if "real_data" in item.keywords:
                item.add_marker(skip_real_data)

    if config.getoption("--no-slow"):
        # Skip test that are marked with `slow` when the command line
        # flag `--no-slow` is provided
        skip_slow = pytest.mark.skip(reason="--no-slow options was given")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)


##################################################################################
# Common Fixtures
##################################################################################
@pytest.fixture(scope="session")
def conda_envs() -> CondaEnv:
    """Returns a ``CondaEnv`` object.

    If a conda environment does not exists then it will build it in a central
    cache. To copy all conda environments into the working directory ``working_dir`` use:

        >>> conda_envs.copy_all_envs(working_dir)
    """
    return CondaEnv()


@pytest.fixture(autouse=True)
def change_default_behavior_of_ConfigMgr(monkeypatch):
    """Monkeypatch the ``ConfigMgr.instance()`` method.

    The ``ConfigMgr`` is designed to create only a single instance. If a
    ``ConfigMgr`` object already exists in the python session it should just
    return that object. This is important to make sure we only have one
    ``ConfigMgr`` object when running the workflow. However, when testing we
    want to create a new object for each test because they all have different
    working directories and configs. Here I monkeypatch the
    ``ConfigMgr.instance()`` method to change it so it always returns a new
    ``ConfigMgr`` object.
    """
    from cgr_gwas_qc.config import ConfigMgr, find_configs

    def mock_instance(*args, **kwargs):
        return ConfigMgr(*find_configs(), **kwargs)

    monkeypatch.setattr(ConfigMgr, "instance", mock_instance)


@pytest.fixture(scope="session")
def qsub():
    """A mock version of qsub.

    Instead of running qsub, just save out the job script and run it with SH.
    """
    return (FakeData() / "scripts/qsub").resolve().as_posix()


##################################################################################
# Small Test Data (Fake)
##################################################################################
@pytest.fixture(scope="session")
def gtc_file() -> Path:
    """Returns the path to a test gtc file."""
    return Path("tests/data/illumina/gtc/small_genotype.gtc").absolute()


@pytest.fixture(scope="session")
def gtc(gtc_file) -> GenotypeCalls:
    """Returns an ``Illumina.GenotypeCalls`` object."""
    return GenotypeCalls(gtc_file)


@pytest.fixture(scope="session")
def vcf_file() -> Path:
    """Returns the path to a test vcf file."""
    return Path("tests/data/1KG/small_1KG.vcf.gz").absolute()


@pytest.fixture(scope="session")
def bpm_file():
    """Returns the path to a test bpm file."""
    return Path("tests/data/illumina/bpm/small_manifest.bpm").absolute()


@pytest.fixture(scope="session")
def bpm(bpm_file) -> BeadPoolManifest:
    """Returns an ``Illumina.BeadPoolManifest`` object."""
    return BeadPoolManifest(bpm_file)


@pytest.fixture(scope="session")
def sample_sheet_file() -> Path:
    """Returns the path to a test sample sheet."""
    return Path("tests/data/example_sample_sheet.csv").absolute()


@pytest.fixture(scope="session")
def sample_sheet(sample_sheet_file) -> SampleSheet:
    """Returns a ``SampleSheet`` object.

    The data section in the sample sheet can be accessed as a
    ``pandas.DataFrame`` using ``sample_sheet.data``.
    """
    return SampleSheet(sample_sheet_file)


@pytest.fixture(scope="session")
def bim_file() -> Path:
    """Returns the path to a test bpm file."""
    return Path("tests/data/plink/samples.bim").absolute()


@pytest.fixture(scope="session")
def vcf(vcf_file) -> pysam.VariantFile:
    """Returns a ``pysam.VariantFile``."""
    return pysam.VariantFile(vcf_file, "r")


@pytest.fixture(scope="session")
def vcf_mock():
    from cgr_gwas_qc.parsers.vcf import VcfFile

    class VcfFileMock(VcfFile):
        def __init__(self, records):
            self.records = records

        def fetch(self, *args, **kwargs):
            return self.records

    return VcfFileMock


##################################################################################
# Update QC Summary Table
##################################################################################
@pytest.mark.real_data
@pytest.fixture(scope="session")
def qc_summary(tmp_path_factory) -> Path:
    """Add new QC columns to legacy QC summary table.

    The legacy QC summary table is missing several columns that are not
    needed by downstream rules. This fixture adds the necessary columns.

    Returns:
        Path to an updated sample qc summary table (CSV).
    """
    from cgr_gwas_qc.workflow.scripts.sample_qc_report import (
        IDENTIFILER_FLAGS,
        QC_HEADER,
        _case_control_encoder,
        _find_study_subject_representative,
        _find_study_subject_with_no_representative,
        _identifiler_reason,
    )

    tmp_path = tmp_path_factory.mktemp("qc_table")
    data_cache = RealData()

    # Add Group By and Internal Control columns
    ss = (
        SampleSheet(data_cache / "original_data/manifest_full.csv")
        .add_group_by_column("PI_Subject_ID")
        .data.assign(Internal_Control=lambda x: x.Sample_Group == "sVALD-001")
        .reindex(
            ["Sample_ID", "Group_By_Subject_ID", "Internal_Control", "Case/Control_Status"], axis=1
        )
    )

    legacy_qc_table = pd.read_csv(data_cache / "production_outputs/all_sample_qc.csv").merge(
        ss, on="Sample_ID", how="left"
    )

    # Use functions from QC script to add other columns
    legacy_qc_table["Case/Control_Status"] = legacy_qc_table["Case/Control_Status"].map(
        _case_control_encoder
    )
    legacy_qc_table["Identifiler_Reason"] = _identifiler_reason(legacy_qc_table, IDENTIFILER_FLAGS)
    legacy_qc_table["Subject_Representative"] = _find_study_subject_representative(legacy_qc_table)
    legacy_qc_table["Subject_Dropped_From_Study"] = _find_study_subject_with_no_representative(
        legacy_qc_table
    )

    legacy_qc_table.reindex(QC_HEADER, axis=1).to_csv(tmp_path / "qc.csv", index=False)
    return tmp_path / "qc.csv"
