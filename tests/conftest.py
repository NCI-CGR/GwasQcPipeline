import shutil
from pathlib import Path

import pytest

from cgr_gwas_qc.config import config_to_yaml
from cgr_gwas_qc.models import config as cfg_models
from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.testing.conda import CondaEnv


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


@pytest.fixture(scope="session")
def conda_envs() -> CondaEnv:
    """Build conda env object.

    If a conda environment does not exists then it will build it. Then I can
    use this object to copy conda environments around my temporary
    test directories.
    """
    return CondaEnv()


def make_config(current_dir: Path, user_files: dict) -> None:
    with chdir(current_dir):
        cfg = cfg_models.Config(
            project_name="Test Project",
            sample_sheet="example_sample_sheet.csv",
            reference_files=cfg_models.ReferenceFiles(
                illumina_manifest_file="small_manifest.bpm",
                thousand_genome_vcf="small_1KG.vcf.gz",
                thousand_genome_tbi="small_1KG.vcf.gz.tbi",
            ),
            user_files=cfg_models.UserFiles(**user_files),
            software_params=cfg_models.SoftwareParams(),
            workflow_params=cfg_models.WorkflowParams(),
        )
        config_to_yaml(cfg)


##################################################################################
# Small Test Data (Fake)
##################################################################################
@pytest.fixture
def working_dir(tmp_path: Path) -> Path:
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
def small_gtc_working_dir(working_dir: Path) -> Path:
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
        shutil.copyfile("tests/data/illumina/gtc/small_genotype.gtc", working_dir / file_name)

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
        shutil.copyfile("tests/data/illumina/idat/small_intensity.idat", working_dir / file_name)

    # create config.yml
    user_files = dict(
        gtc_pattern="{Sample_ID}.gtc",
        idat_pattern=dict(
            red="{SentrixBarcode_A}_{SentrixPosition_A}_Red.idat",
            green="{SentrixBarcode_A}_{SentrixPosition_A}_Grn.idat",
        ),
    )
    make_config(working_dir, user_files)

    return working_dir


@pytest.fixture
def small_ped_working_dir(working_dir: Path) -> Path:
    """Testing working directory with aggregated PED and MAP.

    This is another entry point where the user has pre-aggregated PED and MAP
    files.
    """
    # copy ped/map to working_dir
    shutil.copyfile("tests/data/plink/samples.ped", working_dir / "samples.ped")
    shutil.copyfile("tests/data/plink/samples.map", working_dir / "samples.map")

    # create config.yml
    user_files = dict(ped="samples.ped", map="samples.map")
    make_config(working_dir, user_files)

    return working_dir


@pytest.fixture
def small_bed_working_dir(working_dir: Path) -> Path:
    """Testing working directory with aggregated BED, BIM, and FAM.

    This is another entry point where the user has pre-aggregated BED, BIM,
    and FAM files.
    """
    # copy bed/bim/fam to working_dir
    shutil.copyfile("tests/data/plink/samples.bed", working_dir / "samples.bed")
    shutil.copyfile("tests/data/plink/samples.bim", working_dir / "samples.bim")
    shutil.copyfile("tests/data/plink/samples.fam", working_dir / "samples.fam")

    # create config.yml
    user_files = dict(bed="samples.bed", bim="samples.bim", fam="samples.fam")
    make_config(working_dir, user_files)

    return working_dir
