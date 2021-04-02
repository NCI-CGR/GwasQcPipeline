import importlib.metadata
import os
from pathlib import Path

GITHUB_REF = os.environ.get("GITHUB_REF", None)
if GITHUB_REF and "tags" in GITHUB_REF:
    # On GitHub Actions CI with a tagged release
    __version__ = Path(GITHUB_REF).name  # rev/tags/v0.0.0
else:
    # Use the package version
    __version__ = f"v{importlib.metadata.version('cgr_gwas_qc')}"
