from pathlib import Path
from typing import Optional, Tuple

import pandas as pd
from snakemake.rules import expand

import cgr_gwas_qc.yaml as yaml
from cgr_gwas_qc.models.config import Config
from cgr_gwas_qc.parsers.sample_sheet import SampleSheet


class CgrGwasQcConfigError(Exception):
    pass


class ConfigMgr:
    """Manage all things config.

    There should only be a single instance of this running at a time.

    Attributes:
        SRC_DIR (Path): The absolute path to ``cgr_gwas_qc`` source code.
        WORKFLOW_DIR (Path): The absolute path to the ``workflow`` source code.
        CONFIG_DIR (Path): The absolute path to the ``workflow/config``.
        RULE_DIR (Path): The absolute path to the ``workflow/rules``.
        SCRIPTS_DIR (Path): The absolute path to the ``workflow/scripts``.
        SNAKEFILE (Path): The absolute path to the ``workflow/Snakefile``.

        root (Path): The current working directory.
        user_config (Optional[Path]): The config.yml in the current working directory.
        user_patterns (Optional[Path]): The patterns.yml in the current working directory.

        config (Dict-Like): Workflow configuration settings.
        patterns (Dict-Like): Workflow file name patterns.
    """

    SRC_DIR = Path(__file__).parent.absolute()
    WORKFLOW_DIR = SRC_DIR / "workflow"

    RULE_DIR = WORKFLOW_DIR / "rules"
    SCRIPTS_DIR = WORKFLOW_DIR / "scripts"
    SNAKEFILE = WORKFLOW_DIR / "Snakefile"

    __instance = None

    def __init__(self, root: Path, user_config: Path):
        self.root = root
        self.user_config = user_config

        self._config = Config(**yaml.load(self.user_config))

        self.sample_sheet_file = self.config.reference_paths.sample_sheet
        self._sample_sheet = SampleSheet(self.sample_sheet_file)

    @staticmethod
    def find_configs() -> Tuple[Path, Path]:
        root = Path.cwd().absolute()
        user_config = scan_for_yaml(root, "config")

        if user_config is None:
            raise FileNotFoundError("Please run with a `config.yml` in your working directory.")

        return root, user_config

    @classmethod
    def instance(cls):
        """Returns the active ConfigMgr instance."""
        if cls.__instance is None:
            cls.__instance = cls(*cls.find_configs())
        return cls.__instance

    @property
    def config(self) -> Config:
        return self._config

    @property
    def ss(self) -> pd.DataFrame:
        """Access the sample sheet DataFrame."""
        return self._sample_sheet.data

    def expand(self, pattern, combination=zip):
        return expand(pattern, combination, **self.ss.to_dict())


def scan_for_yaml(base_dir: Path, name: str) -> Optional[Path]:
    """Scans a directory for Yaml configs.

    The Yaml format commonly has ``*.yml`` or ``*.yaml`` file extensions.
    This will search for ``{base_dir}/{name}.{yml,yaml}`` and returns
    the path if present.

    Returns:
        Optiona[Path]: Path to the config file
    """
    if (base_dir / f"{name}.yml").exists():
        return base_dir / f"{name}.yml"

    if (base_dir / f"{name}.yaml").exists():
        return base_dir / f"{name}.yaml"

    return None
