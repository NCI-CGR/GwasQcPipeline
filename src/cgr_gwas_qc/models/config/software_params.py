from pydantic import BaseModel, Field, validator


class SoftwareParams(BaseModel):
    """Software parameters used by various tools in the workflow.

    .. code-block:: yaml

        software_params:
            strand: top

            sample_call_rate_1: 0.8
            snp_call_rate_1: 0.8
            sample_call_rate_2: 0.95
            snp_call_rate_2: 0.95

            intensity_threshold: 6000
            contam_threshold: 0.1
            contam_population: AF

            ld_prune_r2: 0.1
            maf_for_ibd: 0.2
            maf_for_hwe: 0.05

            ibd_pi_hat_min: 0.05
            ibd_pi_hat_max: 1.0
            dup_concordance_cutoff: 0.95
            pi_hat_threshold: 0.2

            autosomal_het_threshold: 0.1

    """

    strand: str = Field(
        "top",
        description="Which Illumina strand to use for genotypes when converting GTC to plink {TOP, FWD, PLUS}. "
        "``workflow/scripts/gtc2plink.py``",
    )

    sample_call_rate_1: float = Field(
        0.8,
        gt=0,
        le=1,
        description="Sample call rate filter 1 threshold. "
        "``plink --mind (1 - sample_call_rate_1``; `reference <https://www.cog-genomics.org/plink/1.9/filter#missing>`__",
    )
    snp_call_rate_1: float = Field(
        0.8,
        gt=0,
        le=1,
        description="SNP call rate filter 1 threshold. "
        "``plink --geno (1 - sample_call_rate_1``; `reference <https://www.cog-genomics.org/plink/1.9/filter#missing>`__",
    )
    sample_call_rate_2: float = Field(
        0.95,
        gt=0,
        le=1,
        description="Sample call rate filter 2 threshold. "
        "``plink --mind (1 - sample_call_rate_1``; `reference <https://www.cog-genomics.org/plink/1.9/filter#missing>`__",
    )
    snp_call_rate_2: float = Field(
        0.95,
        gt=0,
        le=1,
        description="SNP call rate filter 2 threshold. "
        "``plink --geno (1 - sample_call_rate_1``; `reference <https://www.cog-genomics.org/plink/1.9/filter#missing>`__",
    )

    intensity_threshold: int = Field(
        6000,
        gt=0,
        description="Median IDAT intensity threshold used to filter samples during estimating contamination check. "
        "``workflow/scripts/agg_contamination.py``",
    )
    contam_threshold: float = Field(
        0.1,
        gt=0,
        le=1,
        description="%Mix cutoff to consider a sample as contaminated. "
        "``workflow/scripts/agg_contamination.py``",
    )
    contam_population: str = Field(
        "AF",
        description="While population from the 1000 Genomes project to use for B-allele frequencies during contamination testing ."
        "Can be one of {AF, EAS_AF, AMR_AF, AFR_AF, EUR_AF, SAS_AF}. "
        "``workflow/scripts/bpm2abf.py``",
    )

    ld_prune_r2: float = Field(
        0.1,
        gt=0,
        le=1,
        description="The r-squared threshold for LD pruning of SNPS for use with IBD and replicate concordance. "
        "``plink --indep-pairwise 50 5 ld_prune_r2``; `reference <https://www.cog-genomics.org/plink/1.9/ld#indep>`__",
    )
    maf_for_ibd: float = Field(
        0.2,
        gt=0,
        le=1,
        description="The minor allele frequency threshold of SNPS for use with IBD and replicate concordance. "
        "``plink --maf maf_for_ibd``; `reference <https://www.cog-genomics.org/plink/1.9/filter#maf>`__",
    )
    maf_for_hwe: float = Field(
        0.05,
        gt=0,
        le=1,
        description="The minor allele frequency threshold of SNPS for use with population level HWE estimates. "
        "``plink --maf maf_for_hwe``; `reference <https://www.cog-genomics.org/plink/1.9/filter#maf>`__",
    )

    ibd_pi_hat_min: float = Field(
        0.05,
        ge=0,
        lt=1,
        description="The minimum IBD pi hat value to save in the results table. "
        "``plink --genome full --min ibd_pi_hat_min``; `reference <https://www.cog-genomics.org/plink/1.9/ibd>`__",
    )
    ibd_pi_hat_max: float = Field(
        1.0,
        gt=0,
        le=1,
        description="The maximum IBD pi hat value to save in the results table. "
        "``plink --genome full --max ibd_pi_hat_max``; `reference <https://www.cog-genomics.org/plink/1.9/ibd>`__",
    )
    dup_concordance_cutoff: float = Field(
        0.95,
        gt=0,
        le=1,
        description="The concordance threshold to consider two samples as replicates. "
        "``workflow/scripts/concordance_table.py``",
    )
    pi_hat_threshold: float = Field(
        0.2,
        gt=0,
        le=1,
        description="The pi hat threshold to consider two samples as related. "
        "The default of 0.2 reports 1st and 2nd degree relatives. "
        "``workflow/scripts/concordance_table.py``",
    )

    autosomal_het_threshold: float = Field(
        0.1,
        gt=0,
        le=1,
        description="The autosomal heterozygosity F coefficient threshold which to flag subject for removal. "
        "``workflow/scripts/plot_autosomal_heterozygosity.py``",
    )

    @validator("strand")
    def validate_strand(cls, v):
        if v.lower() not in ["top", "fwd", "plus"]:
            raise ValueError("Strand must be one of: TOP, FWD, PLUS")
        return v.lower()

    @staticmethod
    def schema_rst():
        """Tweak schema for rendering in Sphinx."""
        import copy
        import json

        content = copy.deepcopy(SoftwareParams.schema())
        content["title"] = "Software Parameters"

        return json.dumps(content, indent=2)
