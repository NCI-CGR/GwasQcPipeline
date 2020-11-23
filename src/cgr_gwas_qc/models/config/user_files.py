from typing import Optional

from pydantic import BaseModel, Field, FilePath, validator


class UserFiles(BaseModel):
    """File name patterns for user provided data.

    The GWAS QC Pipeline requires GTC or IDAT files.
    """

    idat_pattern: Optional["Idat"]

    gtc_pattern: Optional[str] = Field(
        None,
        description="File name pattern for GTC file. Wildcards are columns in the sample sheet.",
    )

    ped: Optional[FilePath] = Field(
        None, description="Aggregated PED file if sample level GTC files are not available.",
    )
    map: Optional[FilePath] = Field(
        None, description="Aggregated MAP file if sample level GTC files are not available.",
    )

    bed: Optional[FilePath] = Field(
        None, description="Aggregated BED file if sample level GTC files are not available.",
    )
    bim: Optional[FilePath] = Field(
        None, description="Aggregated BIM file if sample level GTC files are not available.",
    )
    fam: Optional[FilePath] = Field(
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
