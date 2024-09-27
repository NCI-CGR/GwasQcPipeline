"""Validate the BGZIP format.

A number of file types use the BGZIP format. Here we focus on the VCF.gz and
VCF.gz.tbi.
"""

from pathlib import Path

from cgr_gwas_qc.exceptions import BgzipMagicNumberError, BgzipTruncatedFileError


def validate(file_name: Path):
    check_magic_number(file_name)
    check_eof_marker(file_name)


def check_magic_number(file_name: Path):
    """Bgzip is a special form of Gzip."""
    with file_name.open("rb") as fh:
        gzip_magic_number = b"\x1f\x8b"
        if fh.read(2) != gzip_magic_number:
            raise BgzipMagicNumberError


def check_eof_marker(file_name: Path):
    """Htslib EOF mark.

    BAM/BCF/VCF.gz/TBI all have a EOF mark.
    """
    with file_name.open("rb") as fh:
        eof_mark = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
        fh.seek(-28, 2)
        if fh.read(28) != eof_mark:
            raise BgzipTruncatedFileError
