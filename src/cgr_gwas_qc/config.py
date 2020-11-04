from pathlib import Path
from typing import Callable, List, Optional, Tuple, Union

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

        config : Workflow configuration settings.
        ss: User's sample sheet data.

    Methods:
        expand: Uses columns from the user's sample sheet to expand a file pattern.
        envmodules: Pulls environmental module information from user's config.
        conda: Creates the full path to conda environment.
        rules: Creates the full path to snakemake rule.
        scripts: Screates the full path to an internal script.
    """

    SRC_DIR = Path(__file__).parent.absolute()
    WORKFLOW_DIR = SRC_DIR / "workflow"

    CONDA_DIR = WORKFLOW_DIR / "conda"
    RULE_DIR = WORKFLOW_DIR / "rules"
    MODULE_DIR = WORKFLOW_DIR / "modules"
    SCRIPTS_DIR = WORKFLOW_DIR / "scripts"
    SNAKEFILE = WORKFLOW_DIR / "Snakefile"

    __instance = None

    ################################################################################
    # Set-up
    ################################################################################
    def __init__(self, root: Path, user_config: Path, validate=True):
        self.root = root
        self.user_config = user_config

        if validate:
            self._config = Config(**yaml.load(self.user_config))
            self.sample_sheet_file = self.config.sample_sheet
            self._sample_sheet = SampleSheet(self.sample_sheet_file)
        else:
            # Force loading without validation. Use this only for debugging
            self._config = Config.construct(**yaml.load(self.user_config))
            try:
                self.sample_sheet_file = self.config.sample_sheet
                self._sample_sheet = SampleSheet(self.sample_sheet_file)
            except (AttributeError, FileNotFoundError):
                pass

    @classmethod
    def instance(cls, validate=True):
        """Returns the active ConfigMgr instance.

        This ensures that only 1 ConfigMgr is created per python session.
        """
        if cls.__instance is None:
            cls.__instance = cls(*find_configs(), validate)
        return cls.__instance

    ################################################################################
    # Access to the user's config and Sample Sheet
    ################################################################################
    @property
    def config(self) -> Config:
        return self._config

    @property
    def ss(self) -> pd.DataFrame:
        """Access the sample sheet DataFrame."""
        return self._sample_sheet.data

    ################################################################################
    # Helper functions for snakemake
    ################################################################################
    def expand(
        self, file_pattern: str, combination: Callable = zip, query: Optional[str] = None
    ) -> List[str]:
        """Use sample sheet columns to fill in file pattern

        Args:
            file_pattern: A snakemake compatable file name pattern.
            combination: The combinatorial method to use (see
              snakemake.rules.expand). Defaults to zip.
            query: A pandas.DataFrame.query to use to filter before
              expansion. Defaults to None.

        Returns:
            A list of expanded file names.
        """
        df = self.ss.query(query) if query else self.ss
        return expand(file_pattern, combination, **df.to_dict("list"))

    def envmodules(self, module_name: str) -> str:
        """Returns the HPC environment module from user's config.

        On HPC systems you can use environment modules instead of conda. The
        specific module versions are stored in the main config.
        """
        modules = self.config.dict().get("env_modules")
        return modules.get(module_name, "") if modules else ""

    def conda(self, file_name: str) -> str:
        """Return path to a conda env file.

        Given a conda env file_name, prepends the full path to that file.
        """
        return (self.CONDA_DIR / file_name).as_posix()

    def rules(self, file_name: str) -> str:
        """Return the path to a rule file.

        Given a rule file_name, prepends the full path to that rule.
        """
        return (self.RULE_DIR / file_name).as_posix()

    def modules(self, file_name: str) -> str:
        """Return the path to a module file.

        Given a rule file_name, prepends the full path to that module.
        """
        return (self.MODULE_DIR / file_name).as_posix()

    def scripts(self, file_name: str) -> str:
        """Return the path to an interal script.

        Given a script file_name, prepends the full path to that script.
        """
        return (self.SCRIPTS_DIR / file_name).as_posix()


################################################################################
# Helper functions for Set-up
################################################################################
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


def find_configs() -> Tuple[Path, Path]:
    root = Path.cwd().absolute()
    user_config = scan_for_yaml(root, "config")

    if user_config is None:
        raise FileNotFoundError("Please run with a `config.yml` in your working directory.")

    return root, user_config


def create_yaml_config(
    project_name: str, sample_sheet: Union[str, Path], yaml_file: str = "config.yml", **kwargs
) -> None:
    cfg = Config(
        project_name=project_name, sample_sheet=sample_sheet, **kwargs
    )
    yaml.write(cfg.to_dict(), yaml_file)
