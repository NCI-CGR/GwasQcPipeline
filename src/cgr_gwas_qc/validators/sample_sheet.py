"""Validate that file is a LIMs Sample Sheet."""
from pathlib import Path

import pandas as pd

from cgr_gwas_qc.exceptions import (
    SampleSheetMissingRequiredColumnsError,
    SampleSheetMissingSectionHeaderError,
    SampleSheetMissingValueRequiredColumnsError,
    SampleSheetNullRowError,
    SampleSheetTruncatedFileError,
)
from cgr_gwas_qc.parsers.sample_sheet import SampleManifest


def validate_manifest(filename: Path, *args):
    data = filename.read_text()
    _check_section_headers(data)
    _check_file_truncation(data)
    df = SampleManifest(filename).data
    _validate(df, *args)
    _check_manifest_null_rows(filename.read_text())


def validate_sample_sheet(filename: Path, *args):
    df = pd.read_csv(filename)
    _validate(df, *args)
    _check_null_rows(filename.read_text())


def _validate(
    df: pd.DataFrame, subject_id_column: str, expected_sex_column: str, case_control_column: str
):
    _check_required_columns(df, subject_id_column, expected_sex_column, case_control_column)
    _check_missing_values_required_columns(df, subject_id_column)


def _check_section_headers(data: str):
    """Checks that the three major section headers are present."""
    missing_headers = []
    if "[Header]" not in data:
        missing_headers.append("Header")

    if "[Manifests]" not in data:
        missing_headers.append("Manifests")

    if "[Data]" not in data:
        missing_headers.append("Data")

    if missing_headers:
        raise SampleSheetMissingSectionHeaderError(missing_headers)


def _check_file_truncation(data: str):
    """Checks for truncated files.

    Makes sure the first and last row have the same number of field
    delimiters (",") and last row as an line termination (\\n)
    """
    rows = data.splitlines()
    if rows[0].count(",") != rows[-1].count(","):
        raise SampleSheetTruncatedFileError


def _check_null_rows(data: str):
    """Checks if there are any empty rows in the sample sheet."""
    rows = data.splitlines()
    num_columns = rows[0].count(",")
    empty_row = "," * num_columns
    if empty_row in rows:
        raise SampleSheetNullRowError


def _check_manifest_null_rows(data: str):
    """Checks if there are any empty rows in the sample sheet."""
    rows = data.splitlines()
    num_columns = rows[0].count(",")
    empty_row = "," * num_columns
    data_section_idx = rows.index("[Data]" + empty_row)
    if empty_row in rows[data_section_idx:]:
        raise SampleSheetNullRowError


def _check_required_columns(df: pd.DataFrame, *args):
    """Checks if essential columns have missing values."""
    required_columns = ["Sample_ID", *args]
    missing_required_columns = [column for column in required_columns if column not in df.columns]

    if missing_required_columns:
        raise SampleSheetMissingRequiredColumnsError(missing_required_columns)


def _check_missing_values_required_columns(df: pd.DataFrame, subject_id_column: str):
    """Checks if essential columns have missing values."""
    essential_columns = ["Sample_ID", subject_id_column]

    col_w_missing_values = [column for column in essential_columns if df[column].isnull().any()]

    if col_w_missing_values:
        raise SampleSheetMissingValueRequiredColumnsError(col_w_missing_values)
