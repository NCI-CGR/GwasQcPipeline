from pathlib import Path

import pytest

from cgr_gwas_qc import load_config
from cgr_gwas_qc.testing import chdir


def test_test_data(test_data, tmp_path):
    (
        test_data.copy_sample_sheet(tmp_path)
        .copy_reference_files(tmp_path)
        .copy_user_files(tmp_path)
        .make_config(tmp_path)
    )

    with chdir(tmp_path):
        cfg = load_config()  # Config working
        assert cfg.ss.Sample_ID.shape[0] == 4  # Reference files working
        assert Path(cfg.config.user_files.bed).exists()  # User files working


def test_test_data_tweak_config(test_data, tmp_path):
    test_data.copy_sample_sheet(tmp_path).make_config(tmp_path, software_params=dict(strand="fwd"))
    with chdir(tmp_path):
        cfg = load_config()  # Config working
        assert cfg.config.software_params.strand == "fwd"


@pytest.mark.real_data
def test_real_data_caching(real_data):
    assert (real_data / "original_data/manifest_full.csv").exists()


@pytest.mark.real_data
def test_copy_file_to_tmp_dir(real_data, tmp_path):
    real_data.copy("original_data/manifest_full.csv", tmp_path / "manifest_full.csv")
    assert (tmp_path / "manifest_full.csv").exists()


@pytest.mark.real_data
def test_copy_folder_to_tmp_dir(real_data, tmp_path):
    real_data.copy("production_outputs/all_contam", tmp_path / "all_contam")
    assert (tmp_path / "all_contam").is_dir() & (tmp_path / "all_contam").exists()
    assert (tmp_path / "all_contam/contam.csv").exists()
