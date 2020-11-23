"""Configuration system data models."""
from logging import getLogger
from typing import Optional

from pydantic import BaseModel, Field, FilePath, validator

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
    project_name: str = Field(..., description="The project title.")
    sample_sheet: FilePath = Field(..., description="Path to the sample manifest from LIMs.")
    reference_files: ReferenceFiles  # Paths to reference files.
    user_files: UserFiles  # Paths to user provided files.
    software_params: SoftwareParams  # Various software parameters.
    workflow_params: WorkflowParams  # Parameters to control how the workflow is run.
    env_modules: Optional[EnvModules] = None  # Use these HPC environmental modules."

    @validator("pipeline_version")
    def validate_pipeline_version(cls, v):
        if v != __version__:
            raise ValueError(
                f"You are running {__version__} of the pipeline, not version {v}. "
                "Either update your config or install a different version of the pipeline."
            )
        return v


Config.update_forward_refs()
