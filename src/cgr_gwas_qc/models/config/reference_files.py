from pathlib import Path
from typing import Optional

import pandas as pd
from pydantic import BaseModel, Field, computed_field, field_validator

from cgr_gwas_qc.exceptions import CsvBpmMissingRequiredColumnsError, CsvBpmMultiGenomeError


class IlluminaCsvBpm(BaseModel):
    """A model for the Illumina CSV Bead Pool Manifest file.

    .. code-block:: yaml

        path: /path/to/csv/file/GSAMD-24v1-0_20011747_A1.csv
        required_columns:
            - GenomeBuild
            - SourceSeq
            - SourceStrand
            - MapInfo
            - Chr
    """

    path: Path
    required_columns: set = {"GenomeBuild", "SourceSeq", "SourceStrand", "MapInfo", "Chr"}

    @computed_field
    def assay_lineno(self) -> int:
        header = (
            pd.read_csv(self.path, nrows=15, usecols=[0], names=["Illumina"]).Illumina == "[Assay]"
        )
        return header[header].index[0] + 1

    @computed_field
    def cols(self) -> list:
        return pd.read_csv(self.path, skiprows=self.assay_lineno, nrows=0).columns.to_list()

    def get_pattern_lineno(self, pattern: str = "[Controls]", max_char_pos: int = 50000) -> int:
        import os

        with open(self.path, "rb") as file:
            file.seek(0, os.SEEK_END)
            total_chars = file.tell()
            position = total_chars
            reverse_line_no = 1
            match = False
            pattern_last_char = pattern[-1].encode("utf-8")
            pattern_length = len(pattern)
            pattern_bytes = pattern.encode("utf-8")
            max_char_pos = max(total_chars - max_char_pos, 0)
            while (position > max_char_pos) and (not match):
                char = file.read(1)
                if char == b"\n":
                    reverse_line_no += 1
                elif char == pattern_last_char:
                    file.seek(-10, os.SEEK_CUR)
                    if pattern_bytes == file.read(pattern_length):
                        match = True
                position -= 1
                file.seek(-2, os.SEEK_CUR)
        if match:
            return reverse_line_no
        else:
            return 0

    @computed_field
    def total_lines(self) -> int:
        with open(self.path, "r") as f:
            return sum(1 for _ in f)

    @computed_field
    def data_lines_count(self) -> int:
        """Count the number of data lines in the CSV BPM file."""
        return self.total_lines - self.assay_lineno - self.get_pattern_lineno() - 1

    @computed_field
    def required_columns_present(
        self,
    ) -> bool:
        """check if the required columns are present in the CSV BPM file."""
        return self.required_columns.issubset(set(self.cols))

    @computed_field
    def genome_builds(self) -> list:
        """Get the list of genome builds from the CSV BPM file."""
        GenomeBuilds = pd.read_csv(
            self.path,
            skiprows=self.assay_lineno,
            nrows=self.data_lines_count,
            usecols=[self.cols.index("GenomeBuild")],
            dtype="category",
        )["GenomeBuild"].cat.categories.to_list()
        return GenomeBuilds

    @computed_field
    def contigs(self) -> list:
        return pd.read_csv(
            self.path, usecols=["Chr"], skiprows=self.assay_lineno, dtype="category", engine="c"
        ).Chr.cat.categories.to_list()

    @computed_field
    def chromosome_x_included(self) -> bool:
        """Check if chromosome X is included in the Illumina CSV BPM."""
        return "X" in self.contigs

    @computed_field
    def chromosome_y_included(self) -> bool:
        """Check if chromosome Y is included in the Illumina CSV BPM."""
        return "Y" in self.contigs

    @field_validator("required_columns_present", mode="after")
    @classmethod
    def check_missing_columns(cls, v, info):
        if v:
            return v
        else:
            cols = set(info.data.get("cols", []))
            missing = cls.required_columns - cols
            raise CsvBpmMissingRequiredColumnsError(missing)

    @field_validator("genome_builds", mode="after")
    @classmethod
    def validate_genome_builds(cls, v, info):
        if len(v) > 1:
            raise CsvBpmMultiGenomeError(v.genome_builds)
        else:
            return v

    @classmethod
    def from_str(cls, path: str) -> "IlluminaCsvBpm":
        return cls(path=Path(path))


class ReferenceFiles(BaseModel):
    """A list of reference files used by the pipeline.

    .. code-block:: yaml

        reference_files:
            illumina_manifest_file: /path/to/bpm/file/GSAMD-24v1-0_20011747_A1.bpm
            illumina_csv_bpm: /path/to/csv/file/GSAMD-24v1-0_20011747_A1.csv
            illumina_cluster_file: /path/to/egt-cluster/file/
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
        description="Path to Reference fasta file to be used to convert gtc to bcf. This could be compressed with bgzip but not gzip.",
    )

    illumina_csv_bpm: Optional["IlluminaCsvBpm"] = Field(
        None,
        description="Path to CSV bead pool manifest provided by Illumina to be used for gtc to bcf conversion. If csv_bpm is not provided, insertions/deletions will be skipped in gtc-to-bcf conversion.",
    )

    @field_validator("illumina_csv_bpm", mode="before")
    @classmethod
    def convert_to_model(cls, v):
        if v is None:
            return None
        elif isinstance(v, str):
            return IlluminaCsvBpm.from_str(v)
        return v

    @staticmethod
    def schema_rst():
        """Tweak schema for rendering in Sphinx."""
        import copy
        import json

        content = copy.deepcopy(ReferenceFiles.schema())
        content["title"] = "Reference Files"

        return json.dumps(content, indent=2)
