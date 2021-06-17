"""Configuration system data models."""
from enum import Enum
from logging import getLogger
from pathlib import Path
from typing import Optional, Sequence

from pydantic import BaseModel, Field, validator

from cgr_gwas_qc.version import __version__

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


class GenomeBuild(str, Enum):
    hg37 = "hg37"
    hg38 = "hg38"


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
    genome_build: GenomeBuild = Field(GenomeBuild.hg37, help="The human genome build.")
    snp_array: Optional[str] = Field(
        None, help="Which SNP array was used. Only used for reporting."
    )
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
