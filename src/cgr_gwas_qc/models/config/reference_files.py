from pydantic import BaseModel, Field, FilePath, validator


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
