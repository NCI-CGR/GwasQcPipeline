import pytest


@pytest.mark.real_data
def test_real_data_caching(real_data_cache):
    assert (real_data_cache / "original_data/manifest_full.csv").exists()


@pytest.mark.real_data
def test_copy_file_to_tmp_dir(real_data_cache, tmp_path):
    real_data_cache.copy("original_data/manifest_full.csv", tmp_path / "manifest_full.csv")
    assert (tmp_path / "manifest_full.csv").exists()


@pytest.mark.real_data
def test_copy_folder_to_tmp_dir(real_data_cache, tmp_path):
    real_data_cache.copy("production_outputs/all_contam", tmp_path / "all_contam")
    assert (tmp_path / "all_contam").is_dir() & (tmp_path / "all_contam").exists()
    assert (tmp_path / "all_contam/contam.csv").exists()
