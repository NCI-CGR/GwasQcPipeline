import shutil
from pathlib import Path

import pytest

from cgr_gwas_qc.config import config_to_yaml
from cgr_gwas_qc.models import config as cfg_models
from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.testing.conda import CondaEnv
from cgr_gwas_qc.testing.data import RealDataCache


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
# Helper Functions
##################################################################################
def make_config(current_dir: Path, config: dict) -> None:
    # adding defaults for software and workflow params
    config["software_params"] = cfg_models.SoftwareParams().dict()
    config["workflow_params"] = cfg_models.WorkflowParams().dict()

    # write out the config to `current_dir/config.yml`
    with chdir(current_dir):
        cfg = cfg_models.Config.parse_obj(config)
        config_to_yaml(cfg)


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


##################################################################################
# Small Test Data (Fake)
##################################################################################
@pytest.fixture(scope="session")
def gtc_file():
    return Path("tests/data/illumina/gtc/small_genotype.gtc")


@pytest.fixture(scope="session")
def vcf_file():
    return Path("tests/data/1KG/small_1KG.vcf.gz")


@pytest.fixture(scope="session")
def bpm_file():
    return Path("tests/data/illumina/bpm/small_manifest.bpm")


@pytest.fixture(scope="session")
def sample_sheet_file():
    return Path("tests/data/example_sample_sheet.csv")


@pytest.fixture(scope="session")
def bim_file():
    return Path("tests/data/plink/samples.bim")


@pytest.fixture(scope="session")
def vcf(vcf_file):
    import pysam

    return pysam.VariantFile(vcf_file, "r")


@pytest.fixture
def small_working_dir(tmp_path: Path) -> Path:
    """Base working directory for fake data.

    Directory small test reference data files.
    """
    # copy reference data
    shutil.copy2("tests/data/example_sample_sheet.csv", tmp_path)
    shutil.copy2("tests/data/illumina/bpm/small_manifest.bpm", tmp_path)
    shutil.copy2("tests/data/1KG/small_1KG.vcf.gz", tmp_path)
    shutil.copy2("tests/data/1KG/small_1KG.vcf.gz.tbi", tmp_path)

    return tmp_path


@pytest.fixture
def small_gtc_working_dir(small_working_dir: Path) -> Path:
    """Testing working directory with per sample GTC files.

    This is the most common situation where the user is start with GTC and
    BPM files.
    """

    # create gtc file for multiple samples
    gtc_samples = [
        "SB00001_PB0001_A01.gtc",
        "SB00002_PB0001_B01.gtc",
        "SB00003_PB0001_C01.gtc",
        "SB00004_PB0001_D01.gtc",
    ]

    for file_name in gtc_samples:
        shutil.copyfile("tests/data/illumina/gtc/small_genotype.gtc", small_working_dir / file_name)

    # create idat files for multiple samples
    idat_samples = [
        "201274900048_R01C01_Red.idat",
        "201274900048_R01C01_Grn.idat",
        "201274900048_R03C01_Red.idat",
        "201274900048_R03C01_Grn.idat",
        "201274900048_R05C01_Red.idat",
        "201274900048_R05C01_Grn.idat",
        "201274900048_R07C01_Red.idat",
        "201274900048_R07C01_Grn.idat",
    ]

    for file_name in idat_samples:
        shutil.copyfile(
            "tests/data/illumina/idat/small_intensity.idat", small_working_dir / file_name
        )

    # create small test data config
    cfg = {
        "project_name": "Small Test Project",
        "sample_sheet": "example_sample_sheet.csv",
        "reference_files": {
            "illumina_manifest_file": "small_manifest.bpm",
            "thousand_genome_vcf": "small_1KG.vcf.gz",
            "thousand_genome_tbi": "small_1KG.vcf.gz.tbi",
        },
        "user_files": {
            "gtc_pattern": "{Sample_ID}.gtc",
            "idat_pattern": {
                "red": "{SentrixBarcode_A}_{SentrixPosition_A}_Red.idat",
                "green": "{SentrixBarcode_A}_{SentrixPosition_A}_Grn.idat",
            },
        },
    }
    make_config(small_working_dir, cfg)

    return small_working_dir


@pytest.fixture
def small_ped_working_dir(small_working_dir: Path) -> Path:
    """Testing working directory with aggregated PED and MAP.

    This is another entry point where the user has pre-aggregated PED and MAP
    files.
    """
    # copy ped/map to small_working_dir
    shutil.copyfile("tests/data/plink/samples.ped", small_working_dir / "samples.ped")
    shutil.copyfile("tests/data/plink/samples.map", small_working_dir / "samples.map")

    # create small test data config
    cfg = {
        "project_name": "Small Test Project",
        "sample_sheet": "example_sample_sheet.csv",
        "reference_files": {
            "illumina_manifest_file": "small_manifest.bpm",
            "thousand_genome_vcf": "small_1KG.vcf.gz",
            "thousand_genome_tbi": "small_1KG.vcf.gz.tbi",
        },
        "user_files": {"ped": "samples.ped", "map": "samples.map"},
    }
    make_config(small_working_dir, cfg)

    return small_working_dir


@pytest.fixture
def small_bed_working_dir(small_working_dir: Path) -> Path:
    """Testing working directory with aggregated BED, BIM, and FAM.

    This is another entry point where the user has pre-aggregated BED, BIM,
    and FAM files.
    """
    # copy bed/bim/fam to small_working_dir
    shutil.copyfile("tests/data/plink/samples.bed", small_working_dir / "samples.bed")
    shutil.copyfile("tests/data/plink/samples.bim", small_working_dir / "samples.bim")
    shutil.copyfile("tests/data/plink/samples.fam", small_working_dir / "samples.fam")

    # create small test data config
    cfg = {
        "project_name": "Small Test Project",
        "sample_sheet": "example_sample_sheet.csv",
        "reference_files": {
            "illumina_manifest_file": "small_manifest.bpm",
            "thousand_genome_vcf": "small_1KG.vcf.gz",
            "thousand_genome_tbi": "small_1KG.vcf.gz.tbi",
        },
        "user_files": {"bed": "samples.bed", "bim": "samples.bim", "fam": "samples.fam"},
    }
    make_config(small_working_dir, cfg)

    return small_working_dir


##################################################################################
# Example Data (From Plink Website)
##################################################################################


##################################################################################
# Real Data (Only use internally)
##################################################################################
@pytest.mark.real_data
@pytest.fixture(scope="session")
def real_data_cache() -> RealDataCache:
    return RealDataCache()
