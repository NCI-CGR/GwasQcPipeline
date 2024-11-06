from pathlib import Path
from typing import Optional

import pandas as pd
from pydantic import BaseModel, Field, validator

from cgr_gwas_qc.exceptions import csvbpmMissingRequiredColumnsError, csvbpmMultiGenomeError


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
        None,
        description="Path to Reference fasta file to be used to convert gtc to bcf. This could be a compressed fasta file and have .bgz or .gz extension.",
    )

    illumina_csv_bpm: Optional[Path] = Field(
        None,
        description="Path to CSV bead pool manifest provided by Illumina to be used for gtc to bcf conversion. If csv_bpm is not provided, insertions/deletions will be skipped in gtc-to-bcf conversion.",
    )

    @validator("illumina_csv_bpm")
    def validate_illumina_csv_bpm(cls, v):

        def get_illumina_csv_assay_lineno(v):
            header = pd.read_csv(v, nrows=15, usecols=[0], names=["Illumina"]).Illumina == "[Assay]"
            return header[header].index[0] + 1

        def get_illumina_csv_cols(v, assay_lineno):
            return pd.read_csv(v, skiprows=assay_lineno, nrows=0).columns

        def get_control_lineno(v):
            import os

            with open(v, "rb") as file:
                file.seek(0, os.SEEK_END)
                total_chars = file.tell()
                position = total_chars
                reverse_line_no = 1
                match = False
                while (position < total_chars & position > total_chars - 50000) & (~match):
                    char = file.read(1)
                    if char == b"\n":
                        reverse_line_no += 1
                    elif char == b"]":
                        file.seek(-10, os.SEEK_CUR)
                        if b"[Controls]" == file.read(10):
                            match = True
                    position -= 1
                    file.seek(-2, os.SEEK_CUR)
            if match:
                return reverse_line_no
            else:
                return 0

        def count_lines(file_path):
            with open(file_path, "r") as f:
                return sum(1 for _ in f)

        if v is None:
            return v

        required_columns = ["GenomeBuild", "SourceSeq", "SourceStrand", "MapInfo", "Chr"]
        assay_lineno = get_illumina_csv_assay_lineno(v)
        columns_in_illumina_csv_bpm = get_illumina_csv_cols(v, assay_lineno)
        missing_columns = set(required_columns) - set(columns_in_illumina_csv_bpm)
        if not len(missing_columns) == 0:
            raise csvbpmMissingRequiredColumnsError(missing_columns)

        GenomeBuilds = pd.read_csv(
            v,
            skiprows=assay_lineno,
            nrows=count_lines(v) - (assay_lineno) - get_control_lineno(v) - 1,
            usecols=[columns_in_illumina_csv_bpm.get_loc("GenomeBuild")],
            dtype="category",
        )["GenomeBuild"].cat.categories.to_list()

        if len(GenomeBuilds) > 1:
            raise csvbpmMultiGenomeError(GenomeBuilds)

        return v

    @staticmethod
    def schema_rst():
        """Tweak schema for rendering in Sphinx."""
        import copy
        import json

        content = copy.deepcopy(ReferenceFiles.schema())
        content["title"] = "Reference Files"

        return json.dumps(content, indent=2)
