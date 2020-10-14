from pathlib import Path

from cgr_gwas_qc.parsers.illumina import GenotypeCalls
from cgr_gwas_qc.validators import GwasQcValidationError


def validate(file_name: Path):
    name = file_name.name

    try:
        # Illumina's parser has a bunch of different error checks, so I am just
        # using those to validate the file. However, I will throw custom errors
        # for clarity.
        GenotypeCalls(file_name.as_posix())
    except Exception as err:
        if err.args[0] == "GTC format error: bad format identifier":
            raise GtcMagicNumberError(name)
        elif err.args[0] == "Unsupported GTC File version":
            raise GtcVersionError(name)
        elif err.args[0] == "GTC file is incomplete":
            raise GtcTuncatedFileError(name)


################################################################################
# Custom Exceptions
################################################################################
class GtcError(GwasQcValidationError):
    pass


class GtcMagicNumberError(GtcError):
    def __init__(self, name):
        message = f"{name} has missing or wrong GTC magic number."
        super().__init__(message)


class GtcVersionError(GtcError):
    def __init__(self, name):
        message = f"{name} has unknown or unsupported BPM version."
        super().__init__(message)


class GtcTuncatedFileError(GtcError):
    def __init__(self, name):
        message = f"{name} is missing the EOF mark."
        super().__init__(message)
