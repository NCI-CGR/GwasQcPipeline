import pytest

from cgr_gwas_qc import load_config
from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.testing.data import FakeData, Graf, RealData


def test_fake_data(tmp_path):
    # GIVEN: working directory with
    #   - sample sheet
    #   - user files (bed entry point)
    #   - reference files (in a different location)
    #   - config
    (
        FakeData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .add_user_files()
        .make_config()
    )

    # WHEN: we load the config in the working dir
    with chdir(tmp_path):
        cfg = load_config()  # Config working

        # THEN:
        # We can access the sample sheet
        assert cfg.ss.Sample_ID.shape[0] == 4

        # User file (BED) is accessible and in the working directory
        assert cfg.config.user_files.bed.exists()
        assert cfg.config.user_files.bed.absolute().parent == tmp_path.absolute()

        # Reference file (BPM) is accessible but not in the working directory
        assert cfg.config.reference_files.illumina_manifest_file.exists()
        assert (
            cfg.config.reference_files.illumina_manifest_file.absolute().parent
            != tmp_path.absolute()
        )


def test_fake_data_tweak_config(tmp_path):
    # GIVEN: working directory with
    #   - sample sheet
    #   - config
    (
        FakeData(tmp_path)
        .add_sample_sheet()
        .make_config(software_params=dict(strand="fwd"))  # WHEN: we change a software setting
    )

    # THEN: that setting is written to the config
    with chdir(tmp_path):
        cfg = load_config()  # Config working
        assert cfg.config.software_params.strand == "fwd"


@pytest.mark.real_data
def test_real_data_caching():
    # GIVEN: a real data repository
    real_data = RealData()

    # THEN: we have access to the files in that repository
    assert (real_data / "original_data/manifest_full.csv").exists()


@pytest.mark.real_data
def test_copy_file_to_tmp_dir(tmp_path):
    # GIVEN: working dir with
    #   - sample sheet
    RealData(tmp_path).add_sample_sheet()

    # THEN: the sample sheet exists
    assert (tmp_path / "sample_sheet.csv").exists()


@pytest.mark.real_data
def test_copy_folder_to_tmp_dir(tmp_path):
    # GIVEN: working dir with
    #   - a specific output folder
    RealData(tmp_path).copy("production_outputs/all_contam", "all_contam")

    # THEN:
    # The output folder is in the working directory
    assert (tmp_path / "all_contam").is_dir() & (tmp_path / "all_contam").exists()
    # The files from that folder were copied too
    assert (tmp_path / "all_contam/contam.csv").exists()


def test_download_graf(tmp_path, monkeypatch):
    """Test GRAF cache mechanism.

    GRAF is not in bioconda, so this is a hack to download GRAF and store it
    in a cache. Here I am testing the download mechanism to make sure I have
    it right.
    """
    # GIVEN: an empty cache folder
    def mock_cache(name):
        pth = tmp_path / name
        pth.mkdir(exist_ok=True, parents=True)
        return pth

    monkeypatch.setattr("cgr_gwas_qc.testing.data._make_cgr_cache", mock_cache)

    # WHEN: I create a Graf instance
    graf = Graf()

    # THEN: the cache should have the graf program downloaded
    assert (graf / "graf").exists() & (graf / "graf").is_file()  # graf executable
    assert (graf / "data").exists() & (graf / "data").is_dir()  # graf data
