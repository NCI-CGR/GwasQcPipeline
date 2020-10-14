from pathlib import Path

from cgr_gwas_qc.parsers.illumina import BeadPoolManifest
from cgr_gwas_qc.validators import GwasQcValidationError


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


################################################################################
# Custom Exceptions
################################################################################


class BpmError(GwasQcValidationError):
    pass


class BpmMagicNumberError(BpmError):
    def __init__(self, name):
        message = f"{name} has missing or wrong BPM magic number."
        super().__init__(message)


class BpmVersionError(BpmError):
    def __init__(self, name):
        message = f"{name} has unknown or unsupported BPM version."
        super().__init__(message)


class BpmNormalizationError(BpmError):
    def __init__(self, name):
        message = f"{name} has an invalid normalization ID."
        super().__init__(message)


class BpmEntryError(BpmError):
    def __init__(self, name):
        message = f"{name} has an invalid number of assay entries."
        super().__init__(message)
