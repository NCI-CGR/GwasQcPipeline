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
    """``config.yml`` data model.

    .. code-block:: yaml

        pipeline_version: v1.0.0
        slurm_partition: defq
        project_name: SR0001-001_1_0000000
        sample_sheet: /path/to/manifest/file/SR0001-001_1_AnalysisManifest_0000000.csv
        genome_build: hg37
        snp_array: GSAMD-24v1-0
        num_samples: 6336
        num_snps: 700078
        reference_files: ...  # Reference file namesapace
        user_files: ...  # User file namespace
        software_params: ...  # Software parameter namespace
        workflow_params: ...  # Workflow parameter namespace
        Sample_IDs_to_remove:

    """

    pipeline_version: str = Field(
        __version__,
        description="The version of the GwasQcPipeline to use. "
        "If you want to use a different version you may just edit this value to match. "
        "However, it is suggested that you re-run the entire pipeline in case there are differences between version.",
    )

    slurm_partition: Optional[str] = Field(
        None, description="Name of the Slurm partition to which jobs will be submitted "
        "when using the --slurm-generic submit option."
    )

    project_name: Optional[str] = Field(
        None, description="The title of the project to use during report generation."
    )

    sample_sheet: Path = Field(
        ...,
        description="The path to the sample sheet (or LIMs manifest). "
        "This is the file referenced during ``cgr config`` and used to generate ``cgr_sample_sheet.csv``.",
    )

    genome_build: GenomeBuild = Field(
        GenomeBuild.hg37,
        description="The human genome build. "
        "This field is not actually used by the workflow. "
        "It is only used during ``cgr config`` to select the VCF file when running on CGEMs/CCAD.",
    )

    snp_array: Optional[str] = Field(
        None, description="Which SNP array was used. Only used for reporting."
    )

    num_samples: int = Field(
        ...,
        description="Number of samples, automatically calculated from the sample sheet during ``cgr config``.",
    )

    num_snps: int = Field(
        ...,
        description="Number of markers, automatically calculated from the ``reference_files.illumina_manifest_file``. "
        "We will attempt to calculate this during ``cgr config`` and ``cgr pre-flight``.",
    )

    reference_files: ReferenceFiles = Field(
        ...,
        description="Reference file namespace. "
        "Reference files include the Illumina provided BPM and the 1000 Genomes VCF.",
    )

    user_files: UserFiles = Field(
        ...,
        description="User file namespace. "
        "User files include the user provided genotype data in IDAT/GTC, PED/MAP, or BED/BIM/FAM formats.",
    )

    software_params: SoftwareParams = Field(
        ...,
        description="Software parameter namespace. "
        "This includes all parameters passed to 3rd party and internal software and scripts.",
    )

    workflow_params: WorkflowParams = Field(
        ...,
        description="Workflow parameter namespace. "
        "This includes all parameters used to control workflow behavior.",
    )

    Sample_IDs_to_remove: Optional[Sequence[str]] = Field(
        None,
        description="A list of Sample_IDs to exclude from QC. "
        "This is the easiest way for a user to exclude specific samples from the GwasQcPipeline. "
        "These samples will be flagged as ``is_user_exclusion`` and will be present in the report, "
        "but they will have not have any results from this analysis.",
    )

    @validator("pipeline_version")
    def validate_pipeline_version(cls, v):
        if v != __version__:
            raise ValueError(
                f"You are running {__version__} of the pipeline, not version {v}. "
                "Either update your config or install a different version of the pipeline."
            )
        return v

    @staticmethod
    def schema_rst():
        """Tweak schema for rendering in Sphinx."""
        import copy
        import json

        content = copy.deepcopy(Config.schema())
        content["title"] = "Top Level Config"

        del content["properties"]["genome_build"]["allOf"]
        content["properties"]["genome_build"]["type"] = "string; [hg37|hg38]"
        del content["definitions"]["GenomeBuild"]

        del content["properties"]["reference_files"]["allOf"]
        content["properties"]["reference_files"]["type"] = "namespace"
        del content["definitions"]["ReferenceFiles"]

        del content["properties"]["user_files"]["allOf"]
        content["properties"]["user_files"]["type"] = "namespace"
        del content["definitions"]["UserFiles"]
        del content["definitions"]["Idat"]

        del content["properties"]["software_params"]["allOf"]
        content["properties"]["software_params"]["type"] = "namespace"
        del content["definitions"]["SoftwareParams"]

        del content["properties"]["workflow_params"]["allOf"]
        content["properties"]["workflow_params"]["type"] = "namespace"
        del content["definitions"]["WorkflowParams"]

        content["properties"]["Sample_IDs_to_remove"]["type"] = "list of strings"
        del content["properties"]["Sample_IDs_to_remove"]["items"]

        del content["definitions"]
        return json.dumps(content, indent=2)


Config.update_forward_refs()
