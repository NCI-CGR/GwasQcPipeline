from pydantic import BaseModel, Field, validator


class SoftwareParams(BaseModel):
    strand: str = Field(
        "TOP",
        description="Which strand to use when converting GTC to plink {TOP, FWD, PLUS}. default TOP",
    )
    samp_cr_1: float = Field(0.8, description="Sample call rate filter 1", gt=0, le=1)
    snp_cr_1: float = Field(0.8, description="Marker call rate filter 1", gt=0, le=1)
    samp_cr_2: float = Field(0.95, description="Sample call rate filter 2", gt=0, le=1)
    snp_cr_2: float = Field(0.95, description="Marker call rate filter 2", gt=0, le=1)

    @validator("strand")
    def validate_strand(cls, v):
        if v.lower() not in ["top", "fwd", "plus"]:
            raise ValueError("Strand must be one of: TOP, FWD, PLUS")
        return v.lower()
