from pathlib import Path

from cgr_gwas_qc.exceptions import (
    BpmEntryError,
    BpmMagicNumberError,
    BpmNormalizationError,
    BpmVersionError,
)
from cgr_gwas_qc.parsers.illumina import BeadPoolManifest


def validate(file_name: Path):
    name = file_name.name

    try:
        # Illumina's parser has a bunch of different error checks, so I am just
        # using those to validate the file. However, I will throw custom errors
        # for clarity.
        BeadPoolManifest(file_name.as_posix())
    except Exception as err:
        if err.args[0].startswith("Invalid BPM format"):
            raise BpmMagicNumberError(name)
        elif err.args[0].startswith("Unknown BPM version") or err.args[0].startswith(
            "Unsupported BPM version"
        ):
            raise BpmVersionError(name)
        elif err.args[0].startswith("Manifest format error: read invalid normalization ID"):
            raise BpmNormalizationError(name)
        elif err.args[0].startswith("Manifest format error: read invalid number of assay entries"):
            raise BpmEntryError(name)
        else:
            raise err
