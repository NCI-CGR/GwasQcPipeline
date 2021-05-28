import pytest

from cgr_gwas_qc.parsers.sample_sheet import SampleManifest
from cgr_gwas_qc.testing.data import FakeData
from cgr_gwas_qc.validators.sample_sheet import (
    SampleSheetMissingRequiredColumnsError,
    SampleSheetMissingSectionHeaderError,
    SampleSheetMissingValueRequiredColumnsError,
    SampleSheetNullRowError,
    SampleSheetTruncatedFileError,
    _check_data_section_truncation,
    _check_manifest_null_rows,
    _check_missing_values_required_columns,
    _check_required_columns,
    _check_section_headers,
    validate_manifest,
)


def test_sample_sheet_good_sample_sheet():
    """Make sure a good sample sheet passes validation."""
    # GIVEN: a good sample sheet
    # WHEN-THEN: we validate it raises no errors
    validate_manifest(
        FakeData._data_path / FakeData._sample_sheet,
        FakeData._subject_id_column,
        FakeData._expected_sex_column,
        FakeData._case_control_column,
    )


################################################################################
# Error if missing section header (MissingSectionHeaderError)
################################################################################

missing_sections = [
    # Missing Header Section
    (
        ",,,\n"
        "[Manifests],,,\n"
        ",,,\n"
        "[Data],,,\n"
        "Sample_ID,col2,col3,col4\n"
        "SP00001,001,002,004\n"
    ),
    # Missing Manifest Section
    (
        "[Header],,,\n"
        ",,,\n"
        ",,,\n"
        "[Data],,,\n"
        "Sample_ID,col2,col3,col4\n"
        "SP00001,001,002,004\n"
    ),
    # Missing Data Section
    ("[Header],,,\n" ",,,\n" "[Manifests],,,\n" ",,,\n"),
]


@pytest.mark.parametrize("data", missing_sections)
def test_missing_sections(data: str):
    # GIVEN: a sample sheet with different sections missing
    # WHEN-THEN: we validate_manifest it raises a missing section error
    with pytest.raises(SampleSheetMissingSectionHeaderError):
        _check_section_headers(data)


################################################################################
# Error if file appears truncated (TruncatedFileError)
################################################################################


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(
            "[Header],,,\n"
            ",,,\n"
            "[Manifests],,,\n"
            ",,,\n"
            "[Data],,,\n"
            "Sample_ID,col2,col3,col4\n"
            "SP00001,002,003",
            id="Missing last field",
        ),
        pytest.param(
            "[Header],,,\n"
            ",,,\n"
            "[Manifests],,,\n"
            ",,,\n"
            "[Data],,,\n"
            "Sample_ID,col2,col3,col4\n"
            "SP00001,002,003,004",
            id="No missing data",
            marks=pytest.mark.xfail(reason="No missing data", strict=True),
        ),
    ],
)
def test_truncated_file(data):
    # GIVEN: a sample sheet with a partial last row (missing \n)
    # WHEN-THEN: we validate_manifest it raises a truncation error
    """If file appears truncated should raise an error."""
    with pytest.raises(SampleSheetTruncatedFileError):
        _check_data_section_truncation(data)


################################################################################
# Error if missing row in data section (NullRowError)
################################################################################

null_data_rows = [
    # Last data row is null
    (
        "[Header],,,\n"
        "[Data],,,\n"
        "Sample_ID,col2,col3,col4\n"
        "SP00001,001,002,004\n"
        "SP00002,001,002,003\n"
        "SP00001,001,002,003\n"
        ",,,\n"
    ),
    # Other data row is null
    (
        "[Header],,,\n"
        "[Data],,,\n"
        "Sample_ID,col2,col3,col4\n"
        "SP00001,001,002,004\n"
        "SP00002,001,002,003\n"
        ",,,\n"
        "SP00001,001,002,003\n"
    ),
]


@pytest.mark.parametrize("data", null_data_rows)
def test_null_row_in_data_section(data):
    """If empty row in the Data section, then raise an error."""
    # GIVEN: a sample sheet with an empty row
    # WHEN-THEN: we validate_manifest it raises a null row error
    with pytest.raises(SampleSheetNullRowError):
        _check_manifest_null_rows(data)


################################################################################
# Do not error if missing row in other sections
################################################################################

null_other_rows = [
    # Null Header row is ok
    (
        "[Header],,,\n"
        ",,,\n"
        "[Data],,,\n"
        "Sample_ID,col2,col3,col4\n"
        "SP00001,001,002,004\n"
        "SP00002,001,002,003\n"
        "SP00001,001,002,003\n"
    ),
    # Null Manifests row is ok
    (
        "[Manifests],,,\n"
        ",,,\n"
        "[Data],,,\n"
        "Sample_ID,col2,col3,col4\n"
        "SP00001,001,002,004\n"
        "SP00002,001,002,003\n"
        "SP00001,001,002,003\n"
    ),
]


@pytest.mark.parametrize("data", null_other_rows)
def test_null_row_in_another_section(data):
    """Empty rows in the other sections should do nothing."""
    # GIVEN: a sample sheet with an empty row in the header of manifests section
    # WHEN-THEN: we validate_manifest it and it raises no errors
    _check_manifest_null_rows(data)


################################################################################
# Error if missing required columns {Sample_ID, LIMs_Individual_ID}
# (MissingRequiredColumnError)
################################################################################

missing_req_column = [
    # Missing Sample_ID
    (
        "[Data],,\n"
        "col1,LIMS_Individual_ID,col3,col4\n"
        "SP00001,001,002,004\n"
        "SP00002,001,002,003\n"
    ),
    # Missing LIMS_Individual_ID
    (
        "[Data],,,\n"
        "Sample_ID,col2,col3,col4\n"
        "SP00001,001,002,004\n"
        "SP00002,,002,003\n"
        "SP00002,001,002,003\n"
    ),
]


@pytest.fixture(params=missing_req_column)
def missing_data_req_column(tmp_path, request):
    """Returns the path to a sample sheet with missing required column.

    Currently, ``Sample_ID`` and ``LIMS_Individual_ID`` are required columns.
    This creates a sample sheet without these columns.
    """
    sample_sheet = tmp_path / "sample_sheet.csv"
    sample_sheet.write_text(request.param)
    return SampleManifest(sample_sheet).data


def test_check_required_columns(missing_data_req_column):
    # GIVEN: a sample sheet with a missing required column
    # WHEN-THEN: we validate_manifest it and it raises a missing column error
    with pytest.raises(SampleSheetMissingRequiredColumnsError):
        _check_required_columns(missing_data_req_column, "LIMS_Individual_ID")


################################################################################
# Error if missing value in required columns {Sample_ID, LIMs_Individual_ID}
# (MissingValueRequiredColumnsError)
################################################################################

missing_values_in_req_column = [
    # Missing value in Sample_ID
    (
        "[Data],,,\n"
        "Sample_ID,LIMS_Individual_ID,col3,col4\n"
        "SP00001,001,002,004\n"
        "SP00002,001,002,003\n"
        ",001,002,003\n"
    ),
]


@pytest.fixture(params=missing_values_in_req_column)
def missing_data_values(tmp_path, request):
    """Returns the path to a sample sheet with missing data in a required column.

    Currently, ``Sample_ID`` and ``LIMS_Individual_ID`` are required columns.
    This creates a sample sheet with missind data in these columns.
    """
    sample_sheet = tmp_path / "sample_sheet.csv"
    sample_sheet.write_text(request.param)
    return SampleManifest(sample_sheet).data


def test_null_columns(missing_data_values):
    # GIVEN: a sample sheet with a missing data in a required column
    # WHEN-THEN: we validate_manifest it and it raises a missing value error
    with pytest.raises(SampleSheetMissingValueRequiredColumnsError):
        _check_missing_values_required_columns(missing_data_values)
