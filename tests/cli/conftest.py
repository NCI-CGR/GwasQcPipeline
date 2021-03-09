import pytest

from cgr_gwas_qc.testing.data import FakeData


@pytest.fixture
def test_data(tmp_path):
    """Returns an test working directory.

    The working directory contains:
        - Sample Sheet
        - Reference Files
            - Illumina Manifest (BPM)
            - Thousand Genome Files (VCF, TBI)
        - User Files
            - Array intensities (IDAT)
            - Genotype Calls (GRC)
        - Config File
    """
    return FakeData(tmp_path).add_sample_sheet().add_user_files(entry_point="gtc").make_config()
