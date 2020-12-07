from cgr_gwas_qc.config import ConfigMgr


def load_config() -> ConfigMgr:
    """Load the config manager."""
    return ConfigMgr.instance()
