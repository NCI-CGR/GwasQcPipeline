"""Configuration system data models."""
from logging import getLogger
from typing import Optional

from pydantic import BaseModel, Field, FilePath, validator

from .env_modules import EnvModules

logger = getLogger(__name__)


class ReferenceFiles(BaseModel):
    """A list of reference files used by the pipeline."""

    illumina_manifest_file: FilePath = Field(
        "/DCEG/CGF/Infinium/Resources/Manifests/GSAMD-Files/build37/GSAMD-24v1-0_20011747_A1.bpm",
        description="Path to the array BPM file.",
    )

    # illumina_cluster_file: FilePath = Field("", description="Path to the array cluster EGT file.")

    thousand_genome_vcf: FilePath = Field(
        "/DCEG/CGF/Bioinformatics/Production/data/thousG/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz",
        description="Path to the 1KG VCF file.",
    )

    thousand_genome_tbi: FilePath = Field(
        "/DCEG/CGF/Bioinformatics/Production/data/thousG/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz.tbi",
        description="Path to the corresponding index for the 1KG VCF file.",
    )

    @validator("illumina_manifest_file")
    def validate_bpm(cls, v):
        assert v.exists()
        assert v.suffix == ".bpm"
        return v

    # @validator("illumina_cluster_filer")
    # def is_egt(cls, v):
    #     return v

    @validator("thousand_genome_vcf")
    def validate_vcf(cls, v):
        assert v.exists()
        assert v.suffix == ".gz"
        return v

    @validator("thousand_genome_tbi")
    def validate_tbi(cls, v):
        assert v.exists()
        assert v.suffix == ".tbi"
        return v


class Idat(BaseModel):
    """Each sample has an Idat file in two color channels."""

    red: str = Field(
        "/DCEG/CGF/Infinium/ScanData/CGF/ByProject/{Project}/{SentrixBarcode_A}/{SentrixBarcode_A}_{SentrixPosition_A}_Red.idat",
        description="File name pattern for Red Channel Idat. Wildcards are columns in the sample sheet.",
    )

    green: str = Field(
        "/DCEG/CGF/Infinium/ScanData/CGF/ByProject/{Project}/{SentrixBarcode_A}/{SentrixBarcode_A}_{SentrixPosition_A}_Grn.idat",
        description="File name pattern for Green Channel Idat. Wildcards are columns in the sample sheet.",
    )

    @validator("*")
    def validate_idat_pattern(cls, v):
        if not v.endswith(".idat"):
            raise ValueError("Idat suffix should be *.idat")

        if "{" not in v or "}" not in v:
            raise ValueError(
                "This dose not look like a file pattern. Add wildcards, corresponding "
                "to column names in the sample sheet, surrounded by '{}'."
            )

        return v


class FilePatterns(BaseModel):
    """File name patterns for user provided data.

    The GWAS QC Pipeline requires GTC or IDAT files.
    """

    gtc: str = Field(
        "/DCEG/CGF/Infinium/ScanData/CGF/ByProject/{Project}/{SentrixBarcode_A}/{SentrixBarcode_A}_{SentrixPosition_A}.gtc",
        description="File name pattern for GTC file. Wildcards are columns in the sample sheet.",
    )

    idat: Idat  # Refers to Idat class above.

    @validator("gtc")
    def validate_gtc_pattern(cls, v):
        if not v.endswith(".gtc"):
            raise ValueError("GTC suffix should be *.gtc")

        if "{" not in v or "}" not in v:
            raise ValueError(
                "This dose not look like a file pattern. Add wildcards, corresponding "
                "to column names in the sample sheet, surrounded by '{}'."
            )

        return v


class Config(BaseModel):
    """Data model for config.yml

    The ``config.yml`` contains all workflow configuration information. This
    is a data model of the config file. It allows for clearly defining a
    schema and implementing data validation.
    """

    project_name: str = Field(..., description="The project title.")
    pipeline_version: str = Field(..., description="The version of the pipeline to use.")
    sample_sheet: FilePath = Field(..., description="Path to the sample manifest from LIMs.")
    reference_files: ReferenceFiles  # Refers to ReferenceFiles above.
    file_patterns: FilePatterns  # Refers to FilePatterns above.
    env_modules: Optional[
        EnvModules
    ]  # Module information if using an HPC with environemntal modules."

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
