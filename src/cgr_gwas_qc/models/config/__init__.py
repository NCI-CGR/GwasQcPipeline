"""Configuration system data models."""
from json import loads
from logging import getLogger
from typing import Optional

from pydantic import BaseModel, Field, FilePath, validator

from cgr_gwas_qc.version import __version__

from .env_modules import EnvModules
from .reference_files import ReferenceFiles
from .software_params import SoftwareParams
from .user_files import UserFiles
from .workflow_params import WorkflowParams

logger = getLogger(__name__)


class Config(BaseModel):
    """Data model for config.yml

    The ``config.yml`` contains all workflow configuration information. This
    is a data model of the config file. It allows for clearly defining a
    schema and implementing data validation.
    """

    project_name: str = Field(..., description="The project title.")
    pipeline_version: str = Field(__version__, description="The version of the pipeline to use.")
    sample_sheet: FilePath = Field(..., description="Path to the sample manifest from LIMs.")
    reference_files: ReferenceFiles = ReferenceFiles()  # Paths to reference files.
    user_files: UserFiles = UserFiles()  # Paths to user provided files.
    software_params: SoftwareParams = SoftwareParams()  # Various software parameters.
    workflow_params: WorkflowParams = WorkflowParams()  # Parameters to control how the workflow is run.
    env_modules: Optional[EnvModules] = EnvModules()  # Use these HPC environmental modules."

    @validator("pipeline_version")
    def validate_pipeline_version(cls, v):
        if v != __version__:
            raise ValueError(
                f"You are running {__version__} of the pipeline, not version {v}. "
                "Either update your config or install a different version of the pipeline."
            )
        return v

    @validator("sample_sheet")
    def validate_sample_sheet(cls, v):
        from cgr_gwas_qc.validators.sample_sheet import SampleSheetNullRowError, validate

        try:
            validate(v)
        except SampleSheetNullRowError:
            # Most of the time I don't care about empty rows. The parser will
            # drop them automatically. Just warn me if there are any.
            logger.warn(f"{v.name} contains empty rows.")

        return v

    def to_dict(self):
        """Create serializable dictionary and remove deprecated options.

        With the .dict() method Path objects were not being serialized to
        strings for writing out. However, the .json() method does serialize
        path objects to strings. This just converts to JSON and then to dict
        for easier export to YAML.
        """

        def remove_deprecated(data):
            cleaned = {}
            for k, v in data.items():
                if isinstance(v, dict):
                    cleaned[k] = remove_deprecated(v)
                elif v != "Deprecated":
                    cleaned[k] = v
            return cleaned

        return remove_deprecated(loads(self.json()))
