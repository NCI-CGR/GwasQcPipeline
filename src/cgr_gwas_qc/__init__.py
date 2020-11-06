from cgr_gwas_qc.config import ConfigMgr


def load_config(validate=True) -> ConfigMgr:
    """Load the config manager."""
    return ConfigMgr.instance(validate)
