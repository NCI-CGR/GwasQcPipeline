from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field


class ReferenceFiles(BaseModel):
    """A list of reference files used by the pipeline.

    .. code-block:: yaml

        reference_files:
            illumina_manifest_file: /path/to/bpm/file/GSAMD-24v1-0_20011747_A1.bpm
            illumina_csv_bpm: /path/to/csv/file/GSAMD-24v1-0_20011747_A1.csv
            thousand_genome_vcf: /path/to/thousand/genome/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz
            thousand_genome_tbi: /path/to/thousand/genome/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz.tbi
            reference_fasta: /path/to/reference/fasta/GCA_000001405.15_GRCh38_full_analysis_set.fna.bgz

    """

    illumina_manifest_file: Optional[Path] = Field(
        None, description="Path to the Illumina provided BPM file."
    )

    illumina_cluster_file: Optional[Path] = Field(
        None, description="Path to the array cluster EGT file."
    )

    thousand_genome_vcf: Path = Field(..., description="Path to the 1000 Genomes VCF file.")

    thousand_genome_tbi: Path = Field(
        ..., description="Path to the corresponding index for the 1000 Genomes VCF file."
    )

    reference_fasta: Optional[Path] = Field(
        None, description="Path to Reference fasta file to be used to convert gtc to bcf"
    )

    illumina_csv_bpm: Optional[Path] = Field(
        None,
        description="Path to CSV bead pool manifest provided by Illumina to be used for gtc to bcf conversion. If csv_bpm is not provided, insertions/deletions will be skipped in gtc-to-bcf conversion.",
    )

    @staticmethod
    def schema_rst():
        """Tweak schema for rendering in Sphinx."""
        import copy
        import json

        content = copy.deepcopy(ReferenceFiles.schema())
        content["title"] = "Reference Files"

        return json.dumps(content, indent=2)
