from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field, root_validator, validator


class UserFiles(BaseModel):
    """A list of user provided files or naming patterns.

    .. code-block:: yaml

        user_files:
            output_pattern: '{prefix}/{file_type}.{ext}'
            idat_pattern:
                red: /expample/pattern/wildcards/are/columns/in/sample_sheet_file/{Project}/{Sample_ID}_Red.idat
                green: /expample/pattern/wildcards/are/columns/in/sample_sheet_file/{Project}/{Sample_ID}_Grn.idat
            gtc_pattern: /expample/pattern/wildcards/are/columns/in/sample_sheet_file/{Project}/{Sample_ID}.gtc

    or

    .. code-block:: yaml

        user_files:
            output_pattern: '{prefix}/{file_type}.{ext}'
            ped: /path/to/samples.ped
            map: /path/to/samples.map

    or

    .. code-block:: yaml

        user_files:
            output_pattern: '{prefix}/{file_type}.{ext}'
            bed: /path/to/samples.bed
            bim: /path/to/samples.bim
            fam: /path/to/samples.fam

    .. note::
        The IDAT/GTC patterns, PED/MAP paths, and BED/BIM/FAM paths are mutually exclusive.
        You should only provide one set of patterns/paths.
    """

    output_pattern: str = Field(
        "{prefix}/{file_type}.{ext}",
        description="File naming pattern for deliverable files. "
        "In general, you should not need to edit this. "
        "However, if you do decide to change this pattern you must keep the ``{prefix}``, ``{file_type}``, and ``{ext}`` patterns. "
        "These are filled in by the workflow, see the delivery sub-workflow for more details.",
    )

    idat_pattern: Optional["Idat"] = Field(
        None,
        description="File naming pattern for IDAT files. "
        "There are two IDAT files for each samples (red and green). "
        "You need to provide a file naming pattern for each IDAT file. "
        "Wildcards are indicated by an ``{}``. "
        "Wildcards have to match column names in the sample sheet exactly.",
    )

    # GTC
    gtc_pattern: Optional[str] = Field(
        None,
        description="File name pattern for GTC file. "
        "Wildcards are indicated by an ``{}``. "
        "Wildcards have to match column names in the sample sheet exactly.",
    )

    # PED/MAP
    ped: Optional[Path] = Field(
        None,
        description="The full path to an aggregated PED file if the sample level GTC files are not available.",
    )
    map: Optional[Path] = Field(
        None,
        description="The full path to an aggregated MAP file if the sample level GTC files are not available.",
    )

    # BED/BIM/FAM
    bed: Optional[Path] = Field(
        None,
        description="The full path to an aggregated BED file if the sample level GTC files are not available.",
    )
    bim: Optional[Path] = Field(
        None,
        description="The full path to an aggregated BIM file if the sample level GTC files are not available.",
    )
    fam: Optional[Path] = Field(
        None,
        description="The full path to an aggregated FAM file if the sample level GTC files are not available.",
    )

    # BCF
    bcf: Optional[Path] = Field(
        None,
        description="The full path to an aggregated BCF/VCF file perferably encoding the GenCall scores.",
    )

    @validator("gtc_pattern")
    def validate_gtc_pattern(cls, v):
        if v is None:
            return v

        if not v.endswith(".gtc"):
            raise ValueError("GTC suffix should be *.gtc")

        if "{" not in v or "}" not in v:
            raise ValueError(
                "This dose not look like a file pattern. Add wildcards, corresponding "
                "to column names in the sample sheet, surrounded by '{}'."
            )

        return v

    @root_validator
    def check_ped_map(cls, values):
        ped, map_ = values.get("ped"), values.get("map")
        if ped is not None and map_ is None:
            raise ValueError("If you provide PED you have to provide a MAP too.")

        if ped is None and map_ is not None:
            raise ValueError("If you provide MAP you have to provide PED too.")

        return values

    @root_validator
    def check_bed_bim_fam(cls, values):
        bed, bim, fam = values.get("bed"), values.get("bim"), values.get("fam")

        if bed is not None and (bim is None or fam is None):
            raise ValueError("If you provide BED you have to provide a BIM and FAM too.")

        if bim is not None and (bed is None or fam is None):
            raise ValueError("If you provide BIM you have to provide a BED and FAM too.")

        if fam is not None and (bed is None or bim is None):
            raise ValueError("If you provide FAM you have to provide a BED and BIM too.")

        return values

    @staticmethod
    def schema_rst():
        """Tweak schema for rendering in Sphinx."""
        import json

        content = UserFiles.schema()
        content["title"] = "User Files"

        del content["properties"]["idat_pattern"]["allOf"]
        content["properties"]["idat_pattern"]["type"] = "namespace"
        del content["definitions"]

        return json.dumps(content, indent=2)


class Idat(BaseModel):
    """Each sample has an Idat file in two color channels."""

    red: str = Field(
        ...,
        description="File name pattern for Red Channel Idat. "
        "Wildcards are indicated by an ``{}``. "
        "Wildcards have to match column names in the sample sheet exactly.",
    )

    green: str = Field(
        ...,
        description="File name pattern for Green Channel Idat. "
        "Wildcards are indicated by an ``{}``. "
        "Wildcards have to match column names in the sample sheet exactly.",
    )

    @validator("*")
    def validate_idat_pattern(cls, v):
        if not v.endswith(".idat"):
            raise ValueError("Idat suffix should be *.idat")

        if "{" not in v or "}" not in v:
            raise ValueError(
                "This dose not look like a file pattern. Add wildcards, corresponding "
                "to column names in the sample sheet, surrounded by '{}'."
            )

        return v


UserFiles.update_forward_refs()
