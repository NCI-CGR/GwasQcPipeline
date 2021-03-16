from pydantic import BaseModel, Field, validator


class SoftwareParams(BaseModel):
    sample_call_rate_1: float = Field(0.8, description="Sample call rate filter 1", gt=0, le=1)
    snp_call_rate_1: float = Field(0.8, description="SNP call rate filter 1", gt=0, le=1)
    sample_call_rate_2: float = Field(0.95, description="Sample call rate filter 2", gt=0, le=1)
    snp_call_rate_2: float = Field(0.95, description="SNP call rate filter 2", gt=0, le=1)

    ld_prune_r2: float = Field(
        0.1,
        description="r-squared threshold for LD pruning of SNPS for use with IBD and concordance.",
        gt=0,
        le=1,
    )
    maf_for_ibd: float = Field(
        0.2,
        description="Minor allele frequency threshold of SNPS for use with IBD and concordance.",
        gt=0,
        le=1,
    )
    maf_for_hwe: float = Field(
        0.05, description="Minor allele frequency threshold of SNPs for use with HWE.", gt=0, le=1
    )
    ibd_pi_hat_min: float = Field(0.05, description="Minimum IBD pi hat threshold.", ge=0, lt=1)
    ibd_pi_hat_max: float = Field(1.0, description="Maximum IBD pi hat threshold.", gt=0, le=1)
    dup_concordance_cutoff: float = Field(
        0.95, description="SNP concordance threshold to call samples as replicates.", gt=0, le=1
    )

    intensity_threshold: int = Field(
        6000, description="IDAT intensity threshold used for estimating contamination.", gt=0
    )
    contam_threshold: float = Field(
        0.1, description="Threshold to call a sample contaminated.", gt=0, le=1
    )
    contam_population: str = Field(
        "AF",
        description="Population from 1KG to use for B allele frequencies during contamination testing {AF, EAS_AF, AMR_AF, AFR_AF, EUR_AF, SAS_AF}.",
    )

    pi_hat_threshold: float = Field(
        0.2,
        description="PI hat threshold used for calling related subjects. The default of 0.2 removes 2nd degree or higher relatives.",
        gt=0,
        le=1,
    )
    autosomal_het_threshold: float = Field(
        0.1,
        description="Autosomal heterozygosity F coefficient threshold to remove subjects.",
        gt=0,
        le=1,
    )

    strand: str = Field(
        "top",
        description="Which strand to use for genotypes when converting GTC to plink {TOP, FWD, PLUS}.",
    )

    @validator("strand")
    def validate_strand(cls, v):
        if v.lower() not in ["top", "fwd", "plus"]:
            raise ValueError("Strand must be one of: TOP, FWD, PLUS")
        return v.lower()
