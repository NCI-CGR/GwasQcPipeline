from pydantic import BaseModel, Field, validator


class SoftwareParams(BaseModel):
    strand: str = Field(
        "TOP",
        description="Which strand to use when converting GTC to plink {TOP, FWD, PLUS}. default TOP",
    )

    @validator("strand")
    def validate_strand(cls, v):
        if v.lower() not in ["top", "fwd", "plus"]:
            raise ValueError("Strand must be one of: TOP, FWD, PLUS")
        return v.lower()
