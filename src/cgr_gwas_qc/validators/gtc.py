from pathlib import Path

from cgr_gwas_qc.exceptions import GtcMagicNumberError, GtcTruncatedFileError, GtcVersionError
from cgr_gwas_qc.parsers.illumina import GenotypeCalls


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
            raise GtcTruncatedFileError(name)
        else:
            raise err
