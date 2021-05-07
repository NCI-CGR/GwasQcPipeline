from pathlib import Path
from typing import Callable, List, Optional, Tuple

import pandas as pd
from snakemake.rules import expand

import cgr_gwas_qc.yaml as yaml
from cgr_gwas_qc.models.config import Config
from cgr_gwas_qc.parsers import sample_sheet


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
        scripts: Creates the full path to an internal script.
    """

    SRC_DIR = Path(__file__).parent.absolute()
    WORKFLOW_DIR = SRC_DIR / "workflow"
    TEMPLATE_DIR = SRC_DIR / "reporting/templates"

    CONDA_DIR = WORKFLOW_DIR / "conda"
    RULE_DIR = WORKFLOW_DIR / "rules"
    MODULE_DIR = WORKFLOW_DIR / "modules"
    SCRIPTS_DIR = WORKFLOW_DIR / "scripts"
    SUBWORKFLOW_DIR = WORKFLOW_DIR / "sub_workflows"
    SNAKEFILE = WORKFLOW_DIR / "Snakefile"

    __instance = None

    ################################################################################
    # Set-up
    ################################################################################
    def __init__(self, root: Path, user_config: Path, sample_sheet_file: Path):
        self.root: Path = root
        self.user_config: Path = user_config
        self.sample_sheet_file: Path = sample_sheet_file

        data = yaml.load(self.user_config)
        self._config = Config.parse_obj(data)
        self._sample_sheet: pd.DataFrame = sample_sheet.read(self.sample_sheet_file)

    @classmethod
    def instance(cls, pytest=False):
        """Returns the active ConfigMgr instance.

        This ensures that only 1 ConfigMgr is created per python session.
        """
        if pytest:
            return cls(*find_cgr_files())

        if cls.__instance is None:
            cls.__instance = cls(*find_cgr_files())

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
        return self._sample_sheet

    ################################################################################
    # Helper functions for snakemake
    ################################################################################
    def expand(
        self, file_pattern: str, combination: Callable = zip, query: Optional[str] = None
    ) -> List[str]:
        """Use sample sheet columns to fill in file pattern

        Args:
            file_pattern: A snakemake compatible file name pattern.
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

    @classmethod
    def conda(cls, filename: str) -> str:
        """Return path to a conda env file.

        Given a conda env file_name, prepends the full path to that file.
        """
        return (cls.CONDA_DIR / filename).as_posix()

    @classmethod
    def rules(cls, filename: str) -> str:
        """Return the path to a rule file.

        Given a rule file_name, prepends the full path to that rule.
        """
        return (cls.RULE_DIR / filename).as_posix()

    @classmethod
    def modules(cls, filename: str) -> str:
        """Return the path to a module file.

        Given a rule file_name, prepends the full path to that module.
        """
        return (cls.MODULE_DIR / filename).as_posix()

    @classmethod
    def scripts(cls, filename: str) -> str:
        """Return the path to an internal script.

        Given a script file_name, prepends the full path to that script.
        """
        return (cls.SCRIPTS_DIR / filename).as_posix()

    @classmethod
    def subworkflow(cls, workflow: str) -> str:
        """Return the path to a sub-workflow.

        Given a sub-workflow name, give the full path to the snakefile.
        """
        return (cls.SUBWORKFLOW_DIR / f"{workflow}.smk").as_posix()

    @property
    def docx_template(self) -> str:
        return (self.TEMPLATE_DIR / "cgr_reference.docx").as_posix()


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


def find_cgr_files() -> Tuple[Path, Path, Path]:
    root = Path.cwd().absolute()
    user_config = scan_for_yaml(root, "config")
    sample_sheet_file = root / "cgr_sample_sheet.csv"

    if user_config is None:
        raise FileNotFoundError("Please run with a `config.yml` in your working directory.")

    if not sample_sheet_file.exists():
        raise FileNotFoundError(
            "Missing `cgr_sample_sheet.csv`. "
            "This file is created during `cgr preflight`, "
            "please run preflight prior to trying to run the workflow."
        )

    return root, user_config, sample_sheet_file


def config_to_yaml(cfg: Config, yaml_file: str = "config.yml", exclude_none=True, **kwargs) -> None:
    def serialize(data):
        """Scan for path objects and convert to strings."""
        cleaned = {}
        for k, v in data.items():
            if isinstance(v, dict):
                cleaned[k] = serialize(v)
            elif isinstance(v, Path):
                cleaned[k] = v.as_posix()
            else:
                cleaned[k] = v
        return cleaned

    yaml.write(serialize(cfg.dict(exclude_none=exclude_none, **kwargs)), yaml_file)
