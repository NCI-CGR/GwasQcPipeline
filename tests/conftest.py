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
def real_tmp_path(tmp_path_factory) -> Path:
    tmp_path = tmp_path_factory.mktemp("real_data")
    (tmp_path / "sample_level").mkdir()
    (tmp_path / "sample_level/ancestry").mkdir()
    (tmp_path / "sample_level/concordance").mkdir()
    (tmp_path / "sample_level/call_rate_2").mkdir()
    (tmp_path / "subject_level").mkdir()
    (tmp_path / "population_level").mkdir()
    (tmp_path / "population_level/European").mkdir()
    return tmp_path


@pytest.fixture(scope="session")
def fake_data_cache():
    return FakeData()


@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_data_cache(pytestconfig):
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    return RealData()


##################################################################################
# Configuration
##################################################################################
@pytest.fixture(scope="session")
def software_params():
    return SoftwareParams()


@pytest.fixture(scope="session")
def workflow_params():
    return WorkflowParams()


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


@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_cfg(tmp_path_factory, pytestconfig) -> ConfigMgr:
    """Real Data config manager object."""
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    tmp_path = tmp_path_factory.mktemp("real_config")
    RealData(tmp_path).make_config().make_cgr_sample_sheet()

    with chdir(tmp_path):
        cfg = load_config(pytest=True)

    return cfg


@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_sample_sheet_csv(real_cfg: ConfigMgr) -> Path:
    return real_cfg.root / "cgr_sample_sheet.csv"


@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_cfg_short(tmp_path_factory, pytestconfig) -> ConfigMgr:
    """Real Data config manager object with IDAT/GTC files."""
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    tmp_path = tmp_path_factory.mktemp("real_config_short")
    RealData(tmp_path, full_sample_sheet=False).add_user_files(
        entry_point="gtc"
    ).make_config().make_cgr_sample_sheet()

    with chdir(tmp_path):
        cfg = load_config(pytest=True)

    return cfg


@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_short_sample_sheet_csv(real_cfg_short: ConfigMgr) -> Path:
    return real_cfg_short.root / "cgr_sample_sheet.csv"


##################################################################################
# New Workflow Outputs
##################################################################################
@pytest.mark.real_data
@pytest.fixture(scope="session")
def sample_concordance_plink(real_cfg: ConfigMgr, real_tmp_path) -> Path:
    """Creates the sample level plink concordance table."""
    from cgr_gwas_qc.workflow.scripts import concordance_table

    # GIVEN: Real lmiss files, a fake mapping of array IDs to 1kg rsIDs, and an
    # outfile.
    outfile = real_tmp_path / "sample_level/concordance/plink.csv"

    data_cache = RealData()
    concordance_table.main(
        data_cache / "production_outputs/ibd/samples.genome",
        real_cfg.config.software_params.dup_concordance_cutoff,
        real_cfg.config.software_params.pi_hat_threshold,
        outfile,
    )

    return outfile


@pytest.mark.real_data
@pytest.fixture(scope="session")
def sample_concordance_graf(sample_concordance_plink, real_tmp_path) -> Path:
    """Forges the sample level graf concordance table.

    I don't want to run GRAF to generate this table, so I am building a fake
    version using PLINK results.
    """
    from cgr_gwas_qc.workflow.scripts import concordance_table

    # GIVEN: PLINK results
    df = concordance_table.read(sample_concordance_plink)

    # WHEN: Simulate GRAF results based on PLINK results.
    outfile = real_tmp_path / "sample_level/concordance/graf.tsv"
    df["HG match"] = 0
    df["HG miss"] = 0
    df["HGMR"] = 0.0
    df["AG match"] = 0
    df["AG miss"] = 0
    df["AGMR"] = 0.0
    df["p_value"] = 0.0
    df["geno relation"] = df.is_ge_concordance.replace({True: "ID", False: "UN"})

    # THEN: Save out the "graf" like table.
    df.rename({"ID1": "sample1", "ID2": "sample2"}, axis=1).drop(
        ["PI_HAT", "concordance", "is_ge_pi_hat", "is_ge_concordance"], axis=1
    ).to_csv(outfile, sep="\t")

    return outfile


@pytest.mark.real_data
@pytest.fixture(scope="session")
def sample_concordance_king(sample_concordance_plink, real_tmp_path) -> Path:
    """Forges the sample level king concordance table.

    I don't want to run KING to generate this table, so I am building a fake
    version using PLINK results.
    """
    from cgr_gwas_qc.workflow.scripts import concordance_table

    # GIVEN: PLINK results
    df = concordance_table.read(sample_concordance_plink)

    # WHEN: Simulate KING results based on PLINK results.
    outfile = real_tmp_path / "sample_level/concordance/king.tsv"
    df["N_SNP"] = 0
    df["HetHet"] = 0
    df["IBS0"] = 0.0
    df["Kinship"] = df.is_ge_concordance.replace({True: 0.4, False: 0.01})

    # THEN: Save out the "king" like table.
    df.drop(["PI_HAT", "concordance", "is_ge_pi_hat", "is_ge_concordance"], axis=1).to_csv(
        outfile, sep="\t"
    )

    return outfile


@pytest.mark.real_data
@pytest.fixture(scope="session")
def sample_concordance_csv(
    real_sample_sheet_csv,
    sample_concordance_plink,
    sample_concordance_graf,
    sample_concordance_king,
    real_tmp_path,
):
    from cgr_gwas_qc.workflow.scripts import sample_concordance

    outfile = real_tmp_path / "sample_level/concordance.csv"
    sample_concordance.main(
        real_sample_sheet_csv,
        sample_concordance_plink,
        sample_concordance_graf,
        sample_concordance_king,
        outfile,
    )
    return outfile


@pytest.mark.real_data
@pytest.fixture(scope="session")
def split_sample_concordance_tables(sample_concordance_csv, real_tmp_path):
    from cgr_gwas_qc.workflow.scripts import split_sample_concordance

    split_sample_concordance.main(
        sample_concordance_csv,
        real_tmp_path / "sample_level/concordance/KnownReplicates.csv",
        real_tmp_path / "sample_level/concordance/InternalQcKnown.csv",
        real_tmp_path / "sample_level/concordance/StudySampleKnown.csv",
        real_tmp_path / "sample_level/concordance/UnknownReplicates.csv",
    )
    return real_tmp_path / "sample_level/concordance"


@pytest.mark.real_data
@pytest.fixture(scope="session")
def snp_qc_csv(real_tmp_path) -> Path:
    """The SNP QC table.

    Return:
        Path to a generated SNP QC table generated from real data. Note: 1Kg
        rsID mapping was faked so all values will be missing.
    """
    from cgr_gwas_qc.workflow.scripts import snp_qc_table

    # GIVEN: Real lmiss files, a fake mapping of array IDs to 1kg rsIDs, and an
    # outfile.
    data_cache = RealData()

    fake_1kg = real_tmp_path / "sample_level/call_rate_2/samples_1kg_rsID.csv"
    fake_1kg.write_text("array_id,thousand_genome_id\n")

    outfile = real_tmp_path / "sample_level/snp_qc.csv"

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
def snp_qc_df(snp_qc_csv) -> pd.DataFrame:
    from cgr_gwas_qc.workflow.scripts import snp_qc_table

    return snp_qc_table.read(snp_qc_csv)


@pytest.mark.real_data
@pytest.fixture(scope="session")
def sample_qc_csv(real_cfg, software_params, real_tmp_path) -> Path:
    """The Sample QC table.

    Return:
        Path to a generated Sample QC table generated from real data. Note:
        we do not have GRAF output from the legacy workflow. Instead, we are
        using the SNPweights output which is automatically handled but will
        give a warning.
    """
    from cgr_gwas_qc.workflow.scripts import sample_qc_table

    data_cache = RealData()
    outfile = real_tmp_path / "sample_level/sample_qc.csv"

    sample_qc_table.main(
        real_cfg.root / "cgr_sample_sheet.csv",
        data_cache / "production_outputs/plink_start/samples_start.imiss",
        data_cache / "production_outputs/plink_filter_call_rate_1/samples_filter1.imiss",
        data_cache / "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss",
        data_cache / "production_outputs/plink_filter_call_rate_1/samples_filter1.sexcheck",
        data_cache / "production_outputs/snpweights/samples.snpweights.csv",
        data_cache / "production_outputs/concordance/KnownReplicates.csv",
        data_cache / "production_outputs/concordance/UnknownReplicates.csv",
        data_cache / "production_outputs/all_contam/contam.csv",
        data_cache / "production_outputs/all_sample_idat_intensity/idat_intensity.csv",
        software_params.dup_concordance_cutoff,
        0.2,
        outfile,
    )

    return outfile


@pytest.mark.real_data
@pytest.fixture
def sample_qc_df(sample_qc_csv) -> pd.DataFrame:
    from cgr_gwas_qc.workflow.scripts import sample_qc_table

    return sample_qc_table.read(sample_qc_csv)


@pytest.mark.real_data
@pytest.fixture(scope="session")
def population_qc_csv(real_cfg: ConfigMgr, real_tmp_path) -> Path:
    """The Population QC table.

    Return:
        Path to a generated Population QC table generated from real data.
    """
    from cgr_gwas_qc.workflow.scripts import population_qc_table

    data_cache = RealData()
    relatives = real_tmp_path / "population_level/European/relatives.csv"
    relatives.write_text("QC_Family_ID,relatives\n")

    outfile = real_tmp_path / "population_level/European/qc.csv"

    population_qc_table.main(
        relatives=relatives,
        pca=data_cache / "production_outputs/pca/EUR_subjects.eigenvec",
        autosomal_het=data_cache
        / "production_outputs/autosomal_heterozygosity/EUR_subjects_qc.het",
        population="EUR",
        threshold=real_cfg.config.software_params.autosomal_het_threshold,
        outfile=outfile,
    )

    return outfile


@pytest.mark.real_data
@pytest.fixture
def population_qc_df(population_qc_csv) -> pd.DataFrame:
    from cgr_gwas_qc.workflow.scripts import population_qc_table

    return population_qc_table.read(population_qc_csv)


@pytest.mark.real_data
@pytest.fixture(scope="session")
def agg_population_qc_csv(sample_qc_csv, population_qc_csv, real_tmp_path) -> Path:
    """Create the aggregated population qc table"""
    from cgr_gwas_qc.workflow.scripts import agg_population_qc_tables

    tables = real_tmp_path / "population_level/per_population_qc.done"
    tables.write_text(population_qc_csv.resolve().as_posix())

    outfile = real_tmp_path / "population_level/population_qc.csv"

    agg_population_qc_tables.main(
        sample_qc_table=sample_qc_csv, population_qc_tables=tables, outfile=outfile,
    )

    return outfile


@pytest.mark.real_data
@pytest.fixture
def agg_population_qc_df(agg_population_qc_csv) -> pd.DataFrame:
    from cgr_gwas_qc.workflow.scripts.agg_population_qc_tables import read_agg_population_qc_tables

    return read_agg_population_qc_tables(agg_population_qc_csv)


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

    outfile = fake_tmp_path / "sample_level/concordance.csv"

    data_cache = FakeData()
    sample_concordance.main(
        fake_sample_sheet_csv,
        data_cache / "cgr/concordance.csv",
        data_cache / "graf/relatedness.tsv",
        data_cache / "king/king.kin0",
        outfile,
    )

    return outfile
