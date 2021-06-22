from enum import Enum
from pathlib import Path
from typing import Callable, List, Optional, Tuple

import pandas as pd
from snakemake.rules import expand

import cgr_gwas_qc.yaml as yaml
from cgr_gwas_qc.models.config import Config
from cgr_gwas_qc.parsers import sample_sheet
from cgr_gwas_qc.typing import PathLike


class CgrGwasQcConfigError(Exception):
    pass


class ConfigMgr:
    """Manage all things config.

    Example:
        >>> from cgr_gwas_qc import load_config
        >>> cfg = load_config()  # must have ``config.yml`` and ``cgr_sample_sheet.csv`` in current directory.

    """

    #: The absolute path to ``cgr_gwas_qc`` source code.
    SRC_DIR: Path = Path(__file__).parent.absolute()

    #: The absolute path to the ``workflow`` source code.
    WORKFLOW_DIR: Path = SRC_DIR / "workflow"

    #: The absolute path to the ``workflow/conda``.
    CONDA_DIR: Path = WORKFLOW_DIR / "conda"

    #: The absolute path to the ``workflow/modules``.
    MODULE_DIR: Path = WORKFLOW_DIR / "modules"

    #: The absolute path to the ``workflow/scripts``.
    SCRIPTS_DIR: Path = WORKFLOW_DIR / "scripts"

    #: The absolute path to the ``workflow/sub_workflows``.
    SUBWORKFLOW_DIR: Path = WORKFLOW_DIR / "sub_workflows"

    #: The absolute path to the ``workflow/Snakefile``.
    SNAKEFILE: Path = WORKFLOW_DIR / "Snakefile"

    #: The absolute path to the ``reporting/templates``.
    TEMPLATE_DIR: Path = SRC_DIR / "reporting/templates"

    __instance = None

    ################################################################################
    # Set-up
    ################################################################################
    def __init__(self, root: Path, user_config: Path, sample_sheet_file: Path):
        #: The current working directory.
        self.root: Path = root

        #: The ``config.yml``` in the current working directory.
        self.user_config: Path = user_config

        #: The ``cgr_sample_sheet.csv``` in the current working directory.
        self.sample_sheet_file: Path = sample_sheet_file

        data = yaml.load(self.user_config)
        self._config = Config.parse_obj(data)
        self._sample_sheet: pd.DataFrame = sample_sheet.read(self.sample_sheet_file)
        self._cluster_groups = sorted(self.ss.cluster_group.unique())

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
        """The user's config settings from ``config.yml``.

        Example:
            >>> from cgr_gwas_qc import load_config
            >>> cfg = load_config()
            >>> cfg.config.software_params.maf_for_ibd
            0.2
        """
        return self._config

    @property
    def ss(self) -> pd.DataFrame:
        """The user's sample sheet data.

        Example:
            >>> cfg = load_config()
            >>> cfg.ss.head(2)
            Sample_ID   Subject_ID    ...
                SP001        SB001    ...
                SP002        SB002    ...
        """
        return self._sample_sheet

    @property
    def cluster_groups(self) -> List[str]:
        """List of the cluster group names from ``cgr_sample_sheet.csv``.

        Example:
            >>> cfg = load_config()
            >>> cfg.cluster_groups
            ["cgroup1", "cgroup2"]
        """
        return self._cluster_groups

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

        Example:
            >>> cfg = load_config()  # assuming ``cgr_sample_sheet.csv`` contains 3 samples SP00{1,2,3}
            >>> cfg.expand("{Sample_ID}.txt")
            ["SP001.txt", "SP002.txt", "SP003.txt"]
            >>> cfg.expand("{Sample_ID}.txt", query="Sample_ID == 'SP002'")
            [SP002.txt"]
        """
        df = self.ss.query(query) if query else self.ss
        return expand(file_pattern, combination, **df.to_dict("list"))

    @classmethod
    def conda(cls, filename: str) -> str:
        """Return path to a conda env file.

        Given a conda env file_name, prepends the full path to that file.
        """
        return (cls.CONDA_DIR / f"{filename}.yml").as_posix()

    @classmethod
    def modules(cls, filename: str) -> str:
        """Return the path to a module file.

        Given a rule file_name, prepends the full path to that module.
        """
        return (cls.MODULE_DIR / f"{filename}.smk").as_posix()

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
        """Return the path to the docx template."""
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


def config_to_yaml(
    cfg: Config, yaml_file: PathLike = "config.yml", exclude_none=True, **kwargs
) -> None:
    def serialize(data):
        """Scan for path objects and convert to strings."""
        cleaned = {}
        for k, v in data.items():
            if isinstance(v, dict):
                cleaned[k] = serialize(v)
            elif isinstance(v, Path):
                cleaned[k] = v.as_posix()
            elif isinstance(v, Enum):
                cleaned[k] = v.value
            else:
                cleaned[k] = v
        return cleaned

    yaml.write(serialize(cfg.dict(exclude_none=exclude_none, **kwargs)), Path(yaml_file))
