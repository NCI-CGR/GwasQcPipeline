from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field


class ReferenceFiles(BaseModel):
    """A list of reference files used by the pipeline."""

    illumina_manifest_file: Optional[Path] = Field(None, description="Path to the array BPM file.")

    illumina_cluster_file: Optional[Path] = Field(
        None, description="Path to the array cluster EGT file."
    )

    thousand_genome_vcf: Optional[Path] = Field(None, description="Path to the 1KG VCF file.")

    thousand_genome_tbi: Optional[Path] = Field(
        None, description="Path to the corresponding index for the 1KG VCF file."
    )
