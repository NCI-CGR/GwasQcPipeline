from pathlib import Path

import pytest


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
