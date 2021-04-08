from cgr_gwas_qc.config import ConfigMgr


def load_config(pytest=False) -> ConfigMgr:
    """Load the config manager."""
    return ConfigMgr.instance(pytest)
