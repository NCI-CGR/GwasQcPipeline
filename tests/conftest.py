from pathlib import Path

import pysam
import pytest

from cgr_gwas_qc.parsers.illumina import BeadPoolManifest, GenotypeCalls
from cgr_gwas_qc.parsers.sample_sheet import SampleSheet
from cgr_gwas_qc.testing.conda import CondaEnv
from cgr_gwas_qc.testing.data import RealDataCache, TestData


##################################################################################
# Pytest Configuration and Hooks
##################################################################################
def pytest_addoption(parser):
    """Add pytest command line options."""
    parser.addoption(
        "--real-data", action="store_true", default=False, help="Run tests with real data."
    )


def pytest_collection_modifyitems(config, items):
    """Modify tests based on command line options.

    Based on an example from the pytest website:
        https://docs.pytest.org/en/6.0.1/example/simple.html#control-skipping-of-tests-according-to-command-line-option
    """
    if not config.getoption("--real-data"):
        # Skip test that are marked with `real_data` unless the command line
        # flag `--real-data` is provided
        skip_real_data = pytest.mark.skip(reason="need --real-data option to run")
        for item in items:
            if "real_data" in item.keywords:
                item.add_marker(skip_real_data)


##################################################################################
# Common Fixtures
##################################################################################
@pytest.fixture(scope="session")
def conda_envs() -> CondaEnv:
    """Build conda env object.

    If a conda environment does not exists then it will build it. Then I can
    use this object to copy conda environments around my temporary
    test directories.
    """
    return CondaEnv()


@pytest.fixture(scope="session")
def test_data() -> TestData:
    return TestData()


@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_data() -> RealDataCache:
    return RealDataCache()


@pytest.fixture(autouse=True)
def mock_ConfigMgr(monkeypatch):
    """Monkeypatch the ``ConfigMgr``.

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


##################################################################################
# Small Test Data (Fake)
##################################################################################
@pytest.fixture(scope="session")
def gtc_file() -> Path:
    return Path("tests/data/illumina/gtc/small_genotype.gtc")


@pytest.fixture(scope="session")
def gtc(gtc_file) -> GenotypeCalls:
    return GenotypeCalls(gtc_file)


@pytest.fixture(scope="session")
def vcf_file() -> Path:
    return Path("tests/data/1KG/small_1KG.vcf.gz")


@pytest.fixture(scope="session")
def bpm_file():
    return Path("tests/data/illumina/bpm/small_manifest.bpm")


@pytest.fixture(scope="session")
def bpm(bpm_file) -> BeadPoolManifest:
    return BeadPoolManifest(bpm_file)


@pytest.fixture(scope="session")
def sample_sheet_file() -> Path:
    return Path("tests/data/example_sample_sheet.csv")


@pytest.fixture(scope="session")
def sample_sheet(sample_sheet_file) -> SampleSheet:
    return SampleSheet(sample_sheet_file)


@pytest.fixture(scope="session")
def bim_file() -> Path:
    return Path("tests/data/plink/samples.bim")


@pytest.fixture(scope="session")
def vcf(vcf_file) -> pysam.VariantFile:
    return pysam.VariantFile(vcf_file, "r")
