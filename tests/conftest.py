import os
from pathlib import Path

import pandas as pd
import pytest

from cgr_gwas_qc import load_config
from cgr_gwas_qc.config import ConfigMgr
from cgr_gwas_qc.models.config import SoftwareParams, WorkflowParams
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


@pytest.fixture(scope="session")
def vcf_mock():
    from cgr_gwas_qc.parsers.vcf import VcfFile

    class VcfFileMock(VcfFile):
        def __init__(self, records):
            self.records = records

        def fetch(self, *args, **kwargs):
            return self.records

    return VcfFileMock


@pytest.fixture
def qsub(monkeypatch):
    """Adds the mock version of qsub/qstat/qacct to the path."""
    qsub_dir = Path("tests/data/scripts").resolve().as_posix()
    monkeypatch.setenv("PATH", qsub_dir, prepend=os.pathsep)


@pytest.fixture(scope="session")
def software_params() -> SoftwareParams:
    return SoftwareParams()


@pytest.fixture(scope="session")
def workflow_params() -> WorkflowParams:
    return WorkflowParams()


##################################################################################
# Fake Data
##################################################################################
@pytest.fixture(scope="session")
def fake_tmp_path(tmp_path_factory) -> Path:
    tmp_path = tmp_path_factory.mktemp("fake_data")
    (tmp_path / "sample_level").mkdir()
    (tmp_path / "sample_level/ancestry").mkdir()
    (tmp_path / "sample_level/concordance").mkdir()
    (tmp_path / "subject_level").mkdir()
    (tmp_path / "population_level").mkdir()
    (tmp_path / "population_level/European").mkdir()
    return tmp_path


@pytest.fixture(scope="session")
def fake_data_cache():
    return FakeData()


@pytest.fixture(scope="session")
def fake_cfg(tmp_path_factory) -> ConfigMgr:
    """Fake Data config manager object."""
    tmp_path = tmp_path_factory.mktemp("fake_config")
    FakeData(tmp_path).add_user_files(entry_point="gtc").make_config().make_cgr_sample_sheet()

    with chdir(tmp_path):
        cfg = load_config(pytest=True)

    return cfg


@pytest.mark.real_data
@pytest.fixture(scope="session")
def fake_sample_sheet_csv(fake_cfg: ConfigMgr) -> Path:
    return fake_cfg.root / "cgr_sample_sheet.csv"


@pytest.fixture(scope="session")
def fake_image(fake_tmp_path):
    from PIL import Image

    outfile = fake_tmp_path / "fake.png"
    img = Image.new("RGB", (1024, 1024), color="gray")
    img.save(outfile)

    return outfile


@pytest.fixture(scope="session")
def fake_sample_concordance_csv(fake_sample_sheet_csv, fake_tmp_path):
    from cgr_gwas_qc.workflow.scripts import sample_concordance

    outfile = fake_tmp_path / "sample_level/concordance/summary.csv"

    data_cache = FakeData()
    sample_concordance.main(
        fake_sample_sheet_csv,
        data_cache / "cgr/concordance.csv",
        data_cache / "graf/relatedness.tsv",
        data_cache / "king/king.kin0",
        outfile,
    )

    return outfile


##################################################################################
# Real Data
##################################################################################
@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_data_cache(pytestconfig) -> RealData:
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    return RealData()


@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_cfg_yml(pytestconfig, real_data_cache) -> Path:
    """Real Data config manager object."""
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    return real_data_cache / "config.yml"


@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_cfg(pytestconfig, real_data_cache) -> ConfigMgr:
    """Real Data config manager object."""
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    with chdir(real_data_cache / "dev_outputs"):
        cfg = load_config(pytest=True)

    return cfg


@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_sample_sheet_csv(pytestconfig, real_data_cache) -> Path:
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    return real_data_cache / "cgr_sample_sheet.csv"


@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_sample_sheet(pytestconfig, real_cfg) -> pd.DataFrame:
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    return real_cfg.ss


@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_cfg_short(pytestconfig, tmp_path_factory) -> ConfigMgr:
    """Real Data config manager object with IDAT/GTC files."""
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    tmp_path = tmp_path_factory.mktemp("real_config_short")
    RealData(tmp_path, full_sample_sheet=False).add_user_files(entry_point="gtc").make_config(
        workflow_params=dict(subject_id_column="PI_Subject_ID"),
        software_params=dict(contam_threshold=0.2),
    ).make_cgr_sample_sheet()

    with chdir(tmp_path):
        cfg = load_config(pytest=True)

    return cfg


@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_short_sample_sheet_csv(pytestconfig, real_cfg_short) -> Path:
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    return real_cfg_short.root / "cgr_sample_sheet.csv"


@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_short_sample_sheet(pytestconfig, real_cfg_short) -> pd.DataFrame:
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    return real_cfg_short.ss


@pytest.mark.real_data
@pytest.fixture
def snp_qc_df(pytestconfig, real_data_cache) -> pd.DataFrame:
    from cgr_gwas_qc.workflow.scripts import snp_qc_table

    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    filename = real_data_cache / "dev_outputs/sample_level/snp_qc.csv"
    return snp_qc_table.read(filename)


@pytest.mark.real_data
@pytest.fixture
def sample_qc_df(pytestconfig, real_data_cache) -> pd.DataFrame:
    from cgr_gwas_qc.workflow.scripts import sample_qc_table

    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    filename = real_data_cache / "dev_outputs/sample_level/sample_qc.csv"
    return sample_qc_table.read(filename)


@pytest.mark.real_data
@pytest.fixture
def subject_qc_df(pytestconfig, real_data_cache) -> pd.DataFrame:
    from cgr_gwas_qc.workflow.scripts import subject_qc_table

    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    filename = real_data_cache / "dev_outputs/subject_level/subject_qc.csv"
    return subject_qc_table.read(filename)


@pytest.mark.real_data
@pytest.fixture
def agg_population_qc_df(pytestconfig, real_data_cache) -> pd.DataFrame:
    from cgr_gwas_qc.workflow.scripts import agg_population_qc_tables

    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    filename = real_data_cache / "dev_outputs/subject_level/population_qc.csv"
    return agg_population_qc_tables.read_agg_population_qc_tables(filename)
