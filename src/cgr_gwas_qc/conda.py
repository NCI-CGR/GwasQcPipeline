import hashlib
import os
from pathlib import Path
from typing import Optional

from cgr_gwas_qc.typing import PathLike

CONDA_EXE = Path(os.getenv("CONDA_EXE", "~/miniconda3/bin/conda")).resolve()
CONDA_BASE_PATH = CONDA_EXE.parents[1]


def conda_activate(env_yaml: PathLike) -> str:
    """Return the conda activation string for use in subprocess."""
    env_path = get_snakemake_conda_env(env_yaml)
    return f"source {CONDA_BASE_PATH / 'etc/profile.d/conda.sh'} && conda activate {env_path}"


def get_snakemake_conda_env(env: PathLike, working_dir: Optional[PathLike] = None) -> Path:
    """Working directory env path.

    For snakemake, the environment path should be the md5sum of the full
    path to `{working_dir}/.snakemake/conda` as well as the content of
    the YAML file. Snakemake uses the first 8 characters of this hash.
    """
    env_path = Path(env).resolve()
    working_dir_path = Path(working_dir).resolve() if working_dir else Path.cwd()
    snakemake_conda_path = working_dir_path / ".snakemake/conda"

    md5 = hashlib.md5()
    md5.update(snakemake_conda_path.as_posix().encode())
    md5.update(env_path.read_bytes())
    full_hash = md5.hexdigest()
    full_hash_path = snakemake_conda_path / full_hash
    short_hash_path = snakemake_conda_path / full_hash[:8]  # for envs pre snakemake 6.0.3
    return short_hash_path if short_hash_path.exists() else full_hash_path
