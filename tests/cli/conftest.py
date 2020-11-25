import pytest


@pytest.fixture(autouse=True)
def mock_ConfigMgr(monkeypatch):
    from cgr_gwas_qc.config import ConfigMgr, find_configs

    def mock_instance(*args, **kwargs):
        return ConfigMgr(*find_configs(), **kwargs)

    monkeypatch.setattr(ConfigMgr, "instance", mock_instance)
