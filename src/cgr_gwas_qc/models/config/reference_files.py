from pydantic import BaseModel, Field, FilePath


class ReferenceFiles(BaseModel):
    """A list of reference files used by the pipeline."""

    illumina_manifest_file: FilePath = Field(..., description="Path to the array BPM file.")

    # illumina_cluster_file: FilePath = Field("", description="Path to the array cluster EGT file.")

    thousand_genome_vcf: FilePath = Field(..., description="Path to the 1KG VCF file.")

    thousand_genome_tbi: FilePath = Field(
        ..., description="Path to the corresponding index for the 1KG VCF file."
    )
