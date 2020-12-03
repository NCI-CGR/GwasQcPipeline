import pytest


@pytest.fixture
def example_working_dir(tmp_path, test_data):
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
    (
        test_data.copy_sample_sheet(tmp_path)
        .copy_reference_files(tmp_path)
        .copy_user_files(tmp_path, entry_point="gtc")
        .make_config(tmp_path)
    )
    return tmp_path
