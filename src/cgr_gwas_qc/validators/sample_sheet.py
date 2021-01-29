"""Validate that file is a LIMs Sample Sheet."""
from pathlib import Path

from cgr_gwas_qc.exceptions import (
    SampleSheetMissingRequiredColumnsError,
    SampleSheetMissingSectionHeaderError,
    SampleSheetMissingValueRequiredColumnsError,
    SampleSheetNullRowError,
    SampleSheetTruncatedFileError,
)
from cgr_gwas_qc.parsers.sample_sheet import SampleSheet

REQUIRED_COLUMNS = ["Sample_ID", "LIMS_Individual_ID"]


def validate(file_name: Path):
    name = file_name.name
    data = file_name.read_text()
    check_section_headers(name, data)
    check_file_truncation(name, data)

    ss = SampleSheet(file_name)
    check_required_columns(name, ss)
    check_missing_values_required_columns(name, ss)

    # Keep this check last, because I usually want to ignore it. This way if
    # all other checks pass I can just catch the exception futher up the call
    # stack.
    check_null_rows(name, data)


def check_section_headers(name: str, data: str):
    """Checks that the three major section headers are present."""
    missing_headers = []
    if "[Header]" not in data:
        missing_headers.append("Header")

    if "[Manifests]" not in data:
        missing_headers.append("Manifests")

    if "[Data]" not in data:
        missing_headers.append("Data")

    if missing_headers:
        raise SampleSheetMissingSectionHeaderError(name, missing_headers)


def check_file_truncation(name: str, data: str):
    """Checks for truncated files.

    Makes sure the first and last row have the same number of field
    delimiters (",") and last row as an line termination (\\n)
    """
    rows = data.splitlines()
    if rows[0].count(",") != rows[-1].count(",") or not data.endswith("\n"):
        raise SampleSheetTruncatedFileError(name)


def check_null_rows(name: str, data: str):
    """Checks if there are any empty rows in the data section."""
    rows = data.splitlines()
    num_delim = rows[0].count(",")
    data_idx = rows.index("[Data]" + "," * num_delim)

    data_rows = rows[data_idx:]
    null_row = "," * num_delim
    if null_row in data_rows:
        raise SampleSheetNullRowError(name)


def check_required_columns(name: str, ss: SampleSheet):
    """Checks if essential columns have missing values."""

    missing_required_columns = [
        column for column in REQUIRED_COLUMNS if column not in ss.data.columns
    ]

    if missing_required_columns:
        raise SampleSheetMissingRequiredColumnsError(name, missing_required_columns)


def check_missing_values_required_columns(name: str, ss: SampleSheet):
    """Checks if essential columns have missing values."""
    essential_columns = ["Sample_ID", "LIMS_Individual_ID"]

    col_w_missing_values = [
        column for column in essential_columns if ss.data[column].isnull().any()
    ]

    if col_w_missing_values:
        raise SampleSheetMissingValueRequiredColumnsError(name, col_w_missing_values)
