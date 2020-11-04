from typing import Optional

from pydantic import BaseModel, Field, validator


class UserFiles(BaseModel):
    """File name patterns for user provided data.

    The GWAS QC Pipeline requires GTC or IDAT files.
    """

    idat_pattern: Optional["Idat"]

    gtc_pattern: Optional[str] = Field(
        "/DCEG/CGF/Infinium/ScanData/CGF/ByProject/{Project}/{SentrixBarcode_A}/{SentrixBarcode_A}_{SentrixPosition_A}.gtc",
        description="File name pattern for GTC file. Wildcards are columns in the sample sheet.",
    )

    ped: Optional[str] = Field(
        None, description="Aggregated PED file if sample level GTC files are not available.",
    )
    map: Optional[str] = Field(
        None, description="Aggregated MAP file if sample level GTC files are not available.",
    )

    bed: Optional[str] = Field(
        None, description="Aggregated BED file if sample level GTC files are not available.",
    )
    bim: Optional[str] = Field(
        None, description="Aggregated BIM file if sample level GTC files are not available.",
    )
    fam: Optional[str] = Field(
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

    @validator("ped")
    def validate_ped(cls, v):
        if v is None:
            return v

        if not v.endswith(".ped"):
            raise ValueError("PED suffix should be *.ped")
        return v

    @validator("map")
    def validate_map(cls, v):
        if v is None:
            return v

        if not v.endswith(".map"):
            raise ValueError("MAP suffix should be *.map")
        return v

    @validator("bed")
    def validate_bed(cls, v):
        if v is None:
            return v

        if not v.endswith(".bed"):
            raise ValueError("BED suffix should be *.bed")
        return v

    @validator("bim")
    def validate_bim(cls, v):
        if v is None:
            return v

        if not v.endswith(".bim"):
            raise ValueError("BIM suffix should be *.bim")
        return v

    @validator("fam")
    def validate_fam(cls, v):
        if v is None:
            return v

        if not v.endswith(".fam"):
            raise ValueError("FAM suffix should be *.fam")
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


UserFiles.update_forward_refs()
