from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field, root_validator, validator


class UserFiles(BaseModel):
    """File name patterns for user provided data.

    The GWAS QC Pipeline requires GTC or IDAT files.
    """

    output_pattern: str = Field(
        "{prefix}/{file_type}.{ext}", description="File naming pattern for deliverable files."
    )

    idat_pattern: Optional["Idat"]

    # GTC
    gtc_pattern: Optional[str] = Field(
        None,
        description="File name pattern for GTC file. Wildcards are columns in the sample sheet.",
    )

    # PED/MAP
    ped: Optional[Path] = Field(
        None, description="Aggregated PED file if sample level GTC files are not available.",
    )
    map: Optional[Path] = Field(
        None, description="Aggregated MAP file if sample level GTC files are not available.",
    )

    # BED/BIM/FAM
    bed: Optional[Path] = Field(
        None, description="Aggregated BED file if sample level GTC files are not available.",
    )
    bim: Optional[Path] = Field(
        None, description="Aggregated BIM file if sample level GTC files are not available.",
    )
    fam: Optional[Path] = Field(
        None, description="Aggregated FAM file if sample level GTC files are not available.",
    )

    @validator("gtc_pattern")
    def validate_gtc_pattern(cls, v):
        if v is None:
            return v

        if not v.endswith(".gtc"):
            raise ValueError("GTC suffix should be *.gtc")

        if "{" not in v or "}" not in v:
            raise ValueError(
                "This dose not look like a file pattern. Add wildcards, corresponding "
                "to column names in the sample sheet, surrounded by '{}'."
            )

        return v

    @root_validator
    def check_ped_map(cls, values):
        ped, map_ = values.get("ped"), values.get("map")
        if ped is not None and map_ is None:
            raise ValueError("If you provide PED you have to provide a MAP too.")

        if ped is None and map_ is not None:
            raise ValueError("If you provide MAP you have to provide PED too.")

        return values

    @root_validator
    def check_bed_bim_fam(cls, values):
        bed, bim, fam = values.get("bed"), values.get("bim"), values.get("fam")

        if bed is not None and (bim is None or fam is None):
            raise ValueError("If you provide BED you have to provide a BIM and FAM too.")

        if bim is not None and (bed is None or fam is None):
            raise ValueError("If you provide BIM you have to provide a BED and FAM too.")

        if fam is not None and (bed is None or bim is None):
            raise ValueError("If you provide FAM you have to provide a BED and BIM too.")

        return values


class Idat(BaseModel):
    """Each sample has an Idat file in two color channels."""

    red: str = Field(
        ...,
        description="File name pattern for Red Channel Idat. Wildcards are columns in the sample sheet.",
    )

    green: str = Field(
        ...,
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


UserFiles.update_forward_refs()
