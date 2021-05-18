"""Configuration system data models."""
from logging import getLogger
from pathlib import Path
from typing import Optional, Sequence

from pydantic import BaseModel, Field, validator

from cgr_gwas_qc.version import __version__

from .env_modules import EnvModules
from .reference_files import ReferenceFiles
from .software_params import SoftwareParams
from .user_files import Idat, UserFiles
from .workflow_params import WorkflowParams

__all__ = [
    "Config",
    "EnvModules",
    "Idat",
    "ReferenceFiles",
    "SoftwareParams",
    "UserFiles",
    "WorkflowParams",
]

logger = getLogger(__name__)


class Config(BaseModel):
    """Data model for config.yml

    The ``config.yml`` contains all workflow configuration information. This
    is a data model of the config file. It allows for clearly defining a
    schema and implementing data validation.
    """

    pipeline_version: str = Field(__version__, description="The version of the pipeline to use.")
    project_name: str = Field("Example Project Name", description="The project title.")
    sample_sheet: Path = Field(
        Path("sample_sheet.csv"), description="Path to the sample manifest from LIMs."
    )
    snp_array: str = Field(..., help="Which SNP array was used. Only used for reporting.")
    num_samples: int = Field(
        ..., help="Number of samples, automatically calculated from the sample sheet."
    )
    num_snps: int = Field(
        ...,
        help="Number of markers, automatically calculated from the `reference_files.illumina_manifest_file",
    )
    reference_files: ReferenceFiles
    user_files: UserFiles = UserFiles()  # Paths to user provided files.
    software_params: SoftwareParams = SoftwareParams()  # Various software parameters.
    workflow_params: WorkflowParams = (
        WorkflowParams()
    )  # Parameters to control how the workflow is run.
    env_modules: Optional[EnvModules]  # Use these HPC environmental modules."
    Sample_IDs_to_remove: Optional[Sequence[str]]

    @validator("pipeline_version")
    def validate_pipeline_version(cls, v):
        if v != __version__:
            raise ValueError(
                f"You are running {__version__} of the pipeline, not version {v}. "
                "Either update your config or install a different version of the pipeline."
            )
        return v


Config.update_forward_refs()
