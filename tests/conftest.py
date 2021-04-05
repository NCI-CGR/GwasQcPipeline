import os
from pathlib import Path

import pandas as pd
import pysam
import pytest

from cgr_gwas_qc import load_config
from cgr_gwas_qc.models.config import Config, Idat, SoftwareParams, WorkflowParams
from cgr_gwas_qc.parsers.illumina import BeadPoolManifest, GenotypeCalls
from cgr_gwas_qc.parsers.sample_sheet import SampleSheet
from cgr_gwas_qc.testing import chdir
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


@pytest.fixture
def qsub(monkeypatch):
    """Adds the mock version of qsub/qstat/qacct to the path."""
    qsub_dir = Path("tests/data/scripts").resolve().as_posix()
    monkeypatch.setenv("PATH", qsub_dir, prepend=os.pathsep)


##################################################################################
# Small Test Data (Fake)
##################################################################################
@pytest.fixture(scope="session")
def idat_file() -> Path:
    """Returns the path to a test idat file."""
    return Path("tests/data/illumina/idat/small_intensity.idat").absolute()


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


@pytest.fixture
def sample_sheet(sample_sheet_file) -> pd.DataFrame:
    """FAke Data sample sheet (4 samples)."""
    return SampleSheet(sample_sheet_file).add_group_by_column().data


@pytest.fixture(scope="session")
def fake_config(tmp_path_factory) -> Config:
    """Fake Data config assuming GTC entrypoint."""
    tmp_path = tmp_path_factory.mktemp("fake_config")

    (FakeData(tmp_path).add_user_files("gtc").make_config())

    with chdir(tmp_path):
        cfg = load_config()

    return cfg.config


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
# Test Data (Real)
##################################################################################
@pytest.mark.real_data
@pytest.fixture
def sample_sheet_short(pytestconfig) -> pd.DataFrame:
    """Real Data short sample sheet (2 samples)."""
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    return SampleSheet(RealData() / "original_data/manifest_short.csv").add_group_by_column().data


@pytest.mark.real_data
@pytest.fixture
def sample_sheet_full(pytestconfig) -> pd.DataFrame:
    """Real Data full sample sheet (203 samples)."""
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    return SampleSheet(RealData() / "original_data/manifest_full.csv").add_group_by_column().data


@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_config(tmp_path_factory, pytestconfig) -> Config:
    """Real Data config assuming GTC entrypoint."""
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    tmp_path = tmp_path_factory.mktemp("real_config")

    (
        RealData(tmp_path)
        .add_user_files("gtc")
        .make_config(workflow_params={"subject_id_to_use": "PI_Subject_ID"})
    )

    with chdir(tmp_path):
        cfg = load_config()

    return cfg.config


##################################################################################
# Update QC Summary Table
##################################################################################
@pytest.mark.real_data
@pytest.fixture(scope="session")
def snp_qc(tmp_path_factory) -> Path:
    """The SNP QC table.

    Return:
        Path to a generated SNP QC table generated from real data. Note: 1Kg
        rsID mapping was faked so all values will be missing.
    """
    from cgr_gwas_qc.workflow.scripts import snp_qc_table

    # GIVEN: Real lmiss files, a fake mapping of array IDs to 1kg rsIDs, and an
    # outfile.
    tmp_path = tmp_path_factory.mktemp("snp_qc")
    data_cache = RealData()

    fake_1kg = tmp_path / "1kg.csv"
    fake_1kg.write_text("array_id,thousand_genome_id\n")

    outfile = tmp_path / "snp_qc.csv"

    # WHEN: I run the main function of the snp qc script.
    snp_qc_table.main(
        data_cache / "production_outputs/plink_start/samples_start.lmiss",
        data_cache / "production_outputs/plink_filter_call_rate_1/samples_filter1.lmiss",
        data_cache / "production_outputs/plink_filter_call_rate_2/samples_filter2.lmiss",
        fake_1kg,
        outfile,
    )

    return outfile


@pytest.mark.real_data
@pytest.fixture
def snp_qc_df(snp_qc) -> pd.DataFrame:
    from cgr_gwas_qc.workflow.scripts.snp_qc_table import read_snp_qc

    return read_snp_qc(snp_qc)


@pytest.mark.real_data
@pytest.fixture(scope="session")
def sample_qc(tmp_path_factory) -> Path:
    """The Sample QC table.

    Return:
        Path to a generated Sample QC table generated from real data. Note:
        we do not have GRAF output from the legacy workflow. Instead, we are
        using the SNPweights output which is automatically handled but will
        give a warning.
    """
    from cgr_gwas_qc.workflow.scripts import sample_qc_table

    tmp_path = tmp_path_factory.mktemp("sample_qc")
    outfile = tmp_path / "sample_qc.csv"

    data_cache = RealData(tmp_path)
    workflow_params = WorkflowParams()
    software_params = SoftwareParams()

    sample_qc_table.main(
        data_cache / "original_data/manifest_full.csv",
        data_cache / "production_outputs/plink_start/samples_start.imiss",
        data_cache / "production_outputs/plink_filter_call_rate_1/samples_filter1.imiss",
        data_cache / "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss",
        data_cache / "production_outputs/plink_filter_call_rate_1/samples_filter1.sexcheck",
        data_cache / "production_outputs/snpweights/samples.snpweights.csv",
        data_cache / "production_outputs/concordance/KnownReplicates.csv",
        data_cache / "production_outputs/concordance/UnknownReplicates.csv",
        data_cache / "production_outputs/all_contam/contam.csv",
        data_cache / "production_outputs/all_sample_idat_intensity/idat_intensity.csv",
        workflow_params.expected_sex_col_name,
        Idat(
            red=data_cache._data_path.as_posix()
            + "/original_data/{SentrixBarcode_A}_{SentrixPosition_A}_Red.idat",
            green=data_cache._data_path.as_posix()
            + "/original_data/{SentrixBarcode_A}_{SentrixPosition_A}_Grn.idat",
        ),
        software_params.dup_concordance_cutoff,
        0.2,
        "PI_Subject_ID",
        list(),
        outfile,
    )

    return outfile


@pytest.mark.real_data
@pytest.fixture
def sample_qc_df(sample_qc) -> pd.DataFrame:
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import read_sample_qc

    return read_sample_qc(sample_qc)


@pytest.mark.real_data
@pytest.fixture(scope="session")
def population_qc(real_config, tmp_path_factory) -> Path:
    """The Population QC table.

    Return:
        Path to a generated Population QC table generated from real data.
    """
    from cgr_gwas_qc.workflow.scripts import population_qc_table

    tmp_path = tmp_path_factory.mktemp("population_qc")
    relatives = tmp_path / "relatives.csv"
    relatives.write_text("QC_Family_ID,relatives\n")

    outfile = tmp_path / "population_qc.csv"

    data_cache = RealData(tmp_path)

    population_qc_table.main(
        relatives=relatives,
        pca=data_cache / "production_outputs/pca/EUR_subjects.eigenvec",
        autosomal_het=data_cache
        / "production_outputs/autosomal_heterozygosity/EUR_subjects_qc.het",
        population="EUR",
        threshold=real_config.software_params.autosomal_het_threshold,
        outfile=outfile,
    )

    return outfile


@pytest.mark.real_data
@pytest.fixture
def population_qc_df(population_qc) -> pd.DataFrame:
    from cgr_gwas_qc.workflow.scripts.population_qc_table import read_population_qc

    return read_population_qc(population_qc)


@pytest.mark.real_data
@pytest.fixture(scope="session")
def agg_population_qc(sample_qc, population_qc, tmp_path_factory) -> Path:
    """Create the aggregated population qc table"""
    from cgr_gwas_qc.workflow.scripts import agg_population_qc_tables

    tmp_path = tmp_path_factory.mktemp("agg_population_qc")
    tables = tmp_path / "population_qc_tables.txt"
    tables.write_text(population_qc.resolve().as_posix())

    outfile = tmp_path / "agg_population_qc.csv"

    agg_population_qc_tables.main(
        sample_qc_table=sample_qc, population_qc_tables=tables, outfile=outfile,
    )

    return outfile


@pytest.mark.real_data
@pytest.fixture
def agg_population_qc_df(agg_population_qc) -> pd.DataFrame:
    from cgr_gwas_qc.workflow.scripts.agg_population_qc_tables import read_agg_population_qc_tables

    return read_agg_population_qc_tables(agg_population_qc)


@pytest.fixture(scope="session")
def fake_image(tmp_path_factory):
    from PIL import Image

    tmp_path = tmp_path_factory.mktemp("images")
    outfile = tmp_path / "fake.png"
    img = Image.new("RGB", (1024, 1024), color="gray")
    img.save(outfile)

    return outfile
