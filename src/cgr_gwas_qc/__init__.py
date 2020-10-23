import importlib.metadata

from cgr_gwas_qc.config import ConfigMgr

__version__ = f"v{importlib.metadata.version('cgr_gwas_qc')}"


def load_config(validate=True) -> ConfigMgr:
    """Load the config manager."""
    return ConfigMgr.instance(validate)
