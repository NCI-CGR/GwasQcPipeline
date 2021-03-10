from pathlib import Path

from cgr_gwas_qc.exceptions import IdatMagicNumberError, IdatTruncatedFileError


def validate(file_name: Path):
    check_magic_number(file_name)
    check_end_file_metadata(file_name)


def check_magic_number(file_name):
    with file_name.open("rb") as fh:
        idat_magic_number = b"IDAT"
        if fh.read(4) != idat_magic_number:
            raise IdatMagicNumberError


def check_end_file_metadata(file_name):
    with file_name.open("rb") as fh:
        fh.seek(-200, 2)
        data = fh.read().decode()
        if "Extract Algorithm" not in data:
            raise IdatTruncatedFileError
