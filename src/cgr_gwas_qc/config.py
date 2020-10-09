from pathlib import Path
from typing import ChainMap, List, Optional, Union

import pandas as pd

import cgr_gwas_qc.yaml as yaml
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

    CONFIG_DIR = WORKFLOW_DIR / "config"
    RULE_DIR = WORKFLOW_DIR / "rules"
    SCRIPTS_DIR = WORKFLOW_DIR / "scripts"
    SNAKEFILE = WORKFLOW_DIR / "Snakefile"

    __instance = None

    def __init__(self, root: Path, user_config: Optional[Path], user_patterns: Optional[Path]):
        self.root = root
        self.user_config = user_config
        self.user_patterns = user_patterns

        self._config = yaml.load([self.user_config, self.CONFIG_DIR / "config.yml"])
        self._patterns = yaml.load([self.user_patterns, self.CONFIG_DIR / "patterns.yml"])

        try:
            self.sample_sheet_file = self.config["sample_sheet"]
            self._sample_sheet = SampleSheet(self.sample_sheet_file)
        except IndexError:
            raise CgrGwasQcConfigError("Please add `sample_sheet` to your config file.")
        except FileNotFoundError:
            raise CgrGwasQcConfigError(f"Sample Sheet File [{self.sample_sheet_file}] not found.")

    @staticmethod
    def find_configs():
        root = Path.cwd().absolute()
        user_config = scan_for_yaml(root, "config")
        user_patterns = scan_for_yaml(root, "patterns")
        return root, user_config, user_patterns

    @classmethod
    def instance(cls):
        """Returns the active ConfigMgr instance."""
        if cls.__instance is None:
            cls.__instance = cls(*cls.find_configs())
        return cls.__instance

    @property
    def config(self) -> ChainMap:
        return self._config

    @property
    def patterns(self) -> ChainMap:
        return self._patterns

    @property
    def ss(self) -> pd.DataFrame:
        """Access the sample sheet DataFrame."""
        return self._sample_sheet.data


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


def flatten_nested(content: Union[list, dict, ChainMap]) -> List[Union[str, int, float]]:
    """Flatten an arbitrary nested dictionary.

    Useful for flattening file patterns from a Yaml config.

    Example:
        >>> content = {"one": ["one_a", "one_b"], "two": "two_a"}
        >>> sorted(flatten_nested(content))
        ["one_a", "one_b", "two_a"]

    """

    def gen(iter: Union[list, dict, ChainMap, str, int, float]):
        if isinstance(iter, (dict, ChainMap)):
            yield from flatten_nested(list(iter.values()))
        elif isinstance(iter, list):
            for item in iter:
                yield from flatten_nested(item)
        else:
            yield iter

    return list(gen(content))
