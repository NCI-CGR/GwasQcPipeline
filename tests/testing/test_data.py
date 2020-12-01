import pytest


@pytest.mark.real_data
def test_real_data_caching(real_data_cache):
    assert (real_data_cache / "original_data/manifest_full.csv").exists()
