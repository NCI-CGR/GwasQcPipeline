from pydantic import BaseModel, Field, validator


class SoftwareParams(BaseModel):
    samp_cr_1: float = Field(0.8, description="Sample call rate filter 1", gt=0, le=1)
    snp_cr_1: float = Field(0.8, description="Marker call rate filter 1", gt=0, le=1)
    samp_cr_2: float = Field(0.95, description="Sample call rate filter 2", gt=0, le=1)
    snp_cr_2: float = Field(0.95, description="Marker call rate filter 2", gt=0, le=1)

    ld_prune_r2: float = Field(0.1, description="r-squred threshold for LD pruning.", gt=0, le=1)
    maf_for_ibd: float = Field(0.2, description="Minor allele frequency threshold", gt=0, le=1)
    ibd_pi_hat_cutoff: float = Field(0.95, description="IBD pi hat threshold.", gt=0, le=1)
    dup_concordance_cutoff: float = Field(
        0.95, description="Duplicate concordance cutoff", gt=0, le=1
    )

    contam_threshold: float = Field(0.2, description="Contamination threshold.", gt=0, le=1)
    contam_population: str = Field(
        "AF", description="Population from 1KG to use for contamination testing."
    )

    pi_hat_theshold: float = Field(0.2, description="pi hat threshold.", gt=0, le=1)
    autosomal_het_theshold: float = Field(
        0.1, description="Autosomal heterozygosity threshold.", gt=0, le=1
    )

    minimum_pop_subjects: int = Field(
        50, description="Minimum number of subjects in a population", gt=0
    )
    control_hwp_threshold: int = Field(50, description="Control samples HWP threshold.", gt=0)

    strand: str = Field(
        "TOP",
        description="Which strand to use when converting GTC to plink {TOP, FWD, PLUS}. default TOP",
    )

    @validator("strand")
    def validate_strand(cls, v):
        if v.lower() not in ["top", "fwd", "plus"]:
            raise ValueError("Strand must be one of: TOP, FWD, PLUS")
        return v.lower()
