import shutil
from pathlib import Path

import pytest

from cgr_gwas_qc.config import create_yaml_config
from cgr_gwas_qc.testing import chdir


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
def gtc_working_dir(tmp_path_factory):
    """Testing working directory with per sample GTC files.

    This is the most common situation where the user is start with GTC and
    BPM files.
    """
    tmp_path = tmp_path_factory.mktemp("data_w_gtc_")

    # copy test data to working dir
    shutil.copy2("tests/data/example_sample_sheet.csv", tmp_path)
    shutil.copy2("tests/data/illumina/bpm/small_manifest.bpm", tmp_path)
    shutil.copy2("tests/data/1KG/small_1KG.vcf.gz", tmp_path)
    shutil.copy2("tests/data/1KG/small_1KG.vcf.gz.tbi", tmp_path)
    shutil.copyfile(
        "tests/data/illumina/gtc/small_genotype.gtc", tmp_path / "SB00001_PB0001_A01.gtc"
    )
    shutil.copyfile(
        "tests/data/illumina/gtc/small_genotype.gtc", tmp_path / "SB00002_PB0001_B01.gtc"
    )
    shutil.copyfile(
        "tests/data/illumina/gtc/small_genotype.gtc", tmp_path / "SB00003_PB0001_C01.gtc"
    )
    shutil.copyfile(
        "tests/data/illumina/gtc/small_genotype.gtc", tmp_path / "SB00004_PB0001_D01.gtc"
    )

    # Create a test config
    with chdir(tmp_path):
        create_yaml_config(
            project_name="Test Project",
            sample_sheet="example_sample_sheet.csv",
            reference_files={
                "illumina_manifest_file": "small_manifest.bpm",
                "thousand_genome_vcf": "small_1KG.vcf.gz",
                "thousand_genome_tbi": "small_1KG.vcf.gz.tbi",
            },
            user_files={
                "gtc_pattern": "{Sample_ID}.gtc",
                "idat_pattern": {"red": "{Sample_ID}_Red.idat", "green": "{Sample_ID}_Grn.idat"},
            },
        )

    return tmp_path
