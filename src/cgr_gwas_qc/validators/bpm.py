from pathlib import Path

from cgr_gwas_qc.exceptions import (
    BpmEntryError,
    BpmMagicNumberError,
    BpmNormalizationError,
    BpmVersionError,
)
from cgr_gwas_qc.parsers.illumina import BeadPoolManifest


def validate(file_name: Path):
    try:
        # Illumina's parser has a bunch of different error checks, so I am just
        # using those to validate the file. However, I will throw custom errors
        # for clarity.
        BeadPoolManifest(file_name.as_posix())
    except Exception as err:
        # Unfortunately, the BeadPoolManifest uses the base Exception. I am
        # splitting them out to custom exceptions for easier handling.

        if len(err.args) > 1 and not isinstance(err.args[0], str):
            # I only expect the Exception to have a message, if not then just
            # re-raise.
            raise err

        if err.args[0].startswith("Invalid BPM format"):
            raise BpmMagicNumberError

        if err.args[0].startswith("Unknown BPM version") or err.args[0].startswith(
            "Unsupported BPM version"
        ):
            raise BpmVersionError

        if err.args[0].startswith("Manifest format error: read invalid normalization ID"):
            raise BpmNormalizationError

        if err.args[0].startswith("Manifest format error: read invalid number of assay entries"):
            raise BpmEntryError

        raise err
