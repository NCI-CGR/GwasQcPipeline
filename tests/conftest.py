from pathlib import Path

import pysam
import pytest

from cgr_gwas_qc.parsers.illumina import BeadPoolManifest, GenotypeCalls
from cgr_gwas_qc.parsers.sample_sheet import SampleSheet
from cgr_gwas_qc.testing.conda import CondaEnv
from cgr_gwas_qc.testing.data import RealData


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
    if config.getoption("--real-data"):
        RealData(sync=True)  # Make sure real data is downloaded
    else:
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
