"""Configuration system data models."""
from pydantic import BaseModel, FilePath, validator


class ReferenceFiles(BaseModel):
    """A number of reference files are required by the pipeline.

    sample_sheet: Valid path to CSV exported from LIMS
    illumina_array_manifest: Valid path to array's BPM
    illumina_cluster_file: Valid path to EGT file
    thousand_genome_vcf: Valid path to 1KG VCF
    thousand_genome_tbi: Valid path to 1KG TBI
    """

    sample_sheet: FilePath
    illumina_array_manifest: FilePath
    # illumina_cluster_file: FilePath
    thousand_genome_vcf: FilePath
    thousand_genome_tbi: FilePath

    @validator("sample_sheet")
    def is_sample_sheet(cls, v):
        data = v.read_text()
        # Containes [Header], [Manifests], [Data]
        if "[Header]" not in data:
            raise ValueError("Missing the 'Header' section")
        if "[Manifests]" not in data:
            raise ValueError("Missing the 'Manifests' section")
        if "[Data]" not in data:
            raise ValueError("Missing the 'Data' section")

        try:
            rows = data.split("\n")
            assert rows[0].count(",") > rows[-1].count(",")
        except AssertionError:
            raise AssertionError("Last row is missing fields.")

        return v

    @validator("illumina_array_manifest")
    def is_bpm(cls, v):
        with v.open("rb") as fh:
            try:
                assert fh.read(3) == b"BPM"
            except AssertionError:
                raise AssertionError("Does not contain BPM magic number.")

        return v

    # @validator("illumina_cluster_filer")
    # def is_egt(cls, v):
    #     return v

    @validator("thousand_genome_vcf", "thousand_genome_tbi")
    def is_bgzip(cls, v):
        with v.open("rb") as fh:
            try:
                gzip_magic_number = b"\x1f\x8b"
                assert fh.read(2) == gzip_magic_number
            except AssertionError:
                raise AssertionError("Missing Gzip magic number.")

            # Note: each record of a VCF/TBI starts with \x42\x43\02 but their
            # location is variable so I have not included their check here for
            # simplicity.

            try:
                eof_magic_number = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
                fh.seek(-28, 2)
                assert fh.read(28) == eof_magic_number
            except AssertionError:
                raise AssertionError("Truncated file, missing EOF magic number.")

        return v


class Idat(BaseModel):
    """Each sample has an Idat file in two color channels."""

    red: str
    green: str

    @validator("*")
    def is_gtc_pattern(cls, v):
        try:
            assert v.endswith(".idat")
        except AssertionError:
            raise AssertionError("IDAT suffix should be *.idat")

        try:
            assert "{" in v and "}" in v
        except AssertionError:
            raise AssertionError(
                "This dose not look like a file pattern. Add wildcards, corresponding "
                "to column names in the sample sheet, surrounded by '{}'."
            )


class UserDataPatterns(BaseModel):
    """File name patterns for user provided data.

    The GWAS QC Pipeline requires GTC or IDAT files.
    """

    gtc: str  # File pattern for *.gtc
    idat: Idat  # File pattern for *.idat

    @validator("gtc")
    def is_gtc_pattern(cls, v):
        try:
            assert v.endswith(".gtc")
        except AssertionError:
            raise AssertionError("GTC suffix should be *.gtc")

        try:
            assert "{" in v and "}" in v
        except AssertionError:
            raise AssertionError(
                "This dose not look like a file pattern. Add wildcards, corresponding "
                "to column names in the sample sheet, surrounded by '{}'."
            )


class Config(BaseModel):
    """Data model for config.yml

    The ``config.yml`` contains all workflow configuration information. This
    is a data model of the config file. It allows for clearly defining a
    schema and implementing data validation.
    """

    project_name: str  # Project title
    reference_paths: ReferenceFiles  # Paths to various reference files
    user_data_patterns: UserDataPatterns  # File name pattern for user data files
