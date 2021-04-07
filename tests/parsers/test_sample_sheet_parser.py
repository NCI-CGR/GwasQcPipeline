"""Testing parsing of the Illumina sample sheet."""
from pathlib import Path
from textwrap import dedent

import pandas as pd
import pytest

from cgr_gwas_qc.parsers.sample_sheet import (
    SampleManifest,
    _convert_to_key_value_pair,
    _strip_terminal_commas,
)

################################################################################
# The sample sheet is a mix of a CSV and INI file. Remove trailing commas for
# cleaner parsing.
################################################################################
strings_with_extra_commas = [
    ("key1,value1,,,,,,,,,,\nkey2,value2,,,,,,,,,,\n", "key1,value1\nkey2,value2\n"),
    ("key1,value1,,,,,,,,,,\n,,,,,,,,,,,\n,,,,,,,,,,,\n", "key1,value1\n"),
]


@pytest.mark.parametrize("data,expected", strings_with_extra_commas)
def test_strip_terminal_commas(data, expected):
    # GIVEN: a multi-lined string with multiple-commas on the end
    # WHEN-THEN: we strip those commas they will match the string without the terminal commas
    assert _strip_terminal_commas(data) == expected


################################################################################
# In the Header and Manifests section we want to convert to key-value pairs.
# This tests the conversion of two comma separated values into a dictionary.
################################################################################
strings_without_extra_commas = [
    ("key1,value1\nkey2,value2\n", {"key1": "value1", "key2": "value2"}),
    ("key1,value1\n", {"key1": "value1"}),
]


@pytest.mark.parametrize("data,expected", strings_without_extra_commas)
def test_convert_to_key_value_pair(data, expected):
    # GIVEN: a multi-lined string with a key-value pair separated by a comma
    # WHEN-THEN: they can be converted to a dictionary of kye-value pairs
    assert _convert_to_key_value_pair(data) == expected


################################################################################
# I try to always provide a pathlib.Path object but make sure it will work with
# when the path is a plain string.
################################################################################
example_sample_sheet = "tests/data/example_sample_sheet.csv"


@pytest.mark.parametrize("file_name", [example_sample_sheet, Path(example_sample_sheet)])
def test_load_str_or_path(file_name):
    """Make sure SampleSheet works with strings or paths."""
    # GIVEN-WHEN: we load an example sample sheet
    sample_sheet = SampleManifest(file_name)

    # THEN: we apture the three section headers.
    assert sorted(sample_sheet._sections.keys()) == sorted(["header", "manifests", "data"])


@pytest.fixture(scope="module")
def sample_sheet_obj(sample_sheet_file) -> SampleManifest:
    """Returns a ``SampleSheet`` object.

    The data section in the sample sheet can be accessed as a
    ``pandas.DataFrame`` using ``sample_sheet.data``.
    """
    return SampleManifest(sample_sheet_file)


################################################################################
# The parsed Header and Manifest should be dictionaries while the Data section
# should be a pandas.DataFrame.
################################################################################
def test_sample_sheet_properties_right_type(sample_sheet_obj: SampleManifest):
    # GIVEN: A parsed sample sheet object
    # THEN: that from each section is the correct type
    # header data is a dictionary of key-value pairs
    assert isinstance(sample_sheet_obj.header, dict)
    # manifests data is a dictionary of key-value pairs
    assert isinstance(sample_sheet_obj.manifests, dict)
    # data data is a dataframe
    assert isinstance(sample_sheet_obj.data, pd.DataFrame)


################################################################################
# Check Header and Manifest sections
################################################################################
def test_header_contains_project_info(sample_sheet_obj: SampleManifest):
    assert "SR0001-001;SR0002-001" == sample_sheet_obj.header["Project Name"]


def test_manifests_contains_bpm_info(sample_sheet_obj: SampleManifest):
    assert "GSAMD-24v1-0" == sample_sheet_obj.manifests["snp_array"]
    assert "GSAMD-24v1-0_20011747_A1.bpm" == sample_sheet_obj.manifests["bpm"]


################################################################################
# Sanity check that Data section gives expected results.
################################################################################
def test_sample_sheet_data(sample_sheet_obj: SampleManifest):
    """Test dataframe functionality."""
    # GIVEN: A parsed sample sheet object
    # THEN: the dataframe from the data section behaves as expected
    # has the same number of rows as the example sample sheet data section
    assert sample_sheet_obj.data.shape[0] == 4
    # Allows querying by different fields and returns the right number of results
    assert sample_sheet_obj.data.query("Identifiler_Sex == 'M'").shape[0] == 2
    assert sample_sheet_obj.data.query("`Case/Control_Status` == 'Case'").shape[0] == 2


################################################################################
# Check that empty rows in the data section are  removed.
################################################################################
def test_empty_row_in_sample_sheet_data(tmp_path):
    # GIVEN: A data section with an empty row
    sample_sheet = Path(tmp_path) / "sample_sheet.csv"
    sample_sheet.write_text(
        dedent(
            """\
        [Data],,,
        Sample_ID,col2,col3,col4
        SB001,001,002,004
        SB002,001,002,003
        ,,,
        """
        )
    )

    # WHEN: we parse that sample sheet
    ss = SampleManifest(sample_sheet)

    # THEN: we automatically drop that empty row
    assert ss.data.shape[0] == 2


def test_add_group_by_column(tmp_path):
    # GIVEN: A LIMS sample sheet with the `Group_By` column
    (tmp_path / "sample_sheet.csv").write_text(
        dedent(
            """\
            [Header],,
            test,data,
            [Manifests],,
            test,data,
            [Data]
            Sample_ID,LIMS_Individual_ID,Group_By
            T0001,L0001,Sample_ID
            T0002,L0002,LIMS_Individual_ID
            """
        )
    )
    # WHEN: I parse the sample sheet and add the grouping column
    sample_sheet = SampleManifest(tmp_path / "sample_sheet.csv").add_group_by_column()

    # THEN: A new column `Group_By_Subject_ID` will have the values from the
    # columns specified in the `Group_By` column.
    assert all(sample_sheet.data.Group_By_Subject_ID == ["T0001", "L0002"])


def test_add_user_provided_group_by_column(tmp_path):
    # GIVEN: A LIMS sample sheet with the `Group_By` column
    (tmp_path / "sample_sheet.csv").write_text(
        dedent(
            """\
            [Header],,
            test,data,
            [Manifests],,
            test,data,
            [Data]
            Sample_ID,LIMS_Individual_ID,Group_By
            T0001,L0001,Sample_ID
            T0002,L0002,LIMS_Individual_ID
            """
        )
    )
    # WHEN: I parse the sample sheet and add the grouping column with a user selected column.
    sample_sheet = SampleManifest(tmp_path / "sample_sheet.csv").add_group_by_column("Sample_ID")

    # THEN: A new column `Group_By_Subject_ID` will have the values from the
    # columns specified by the user.
    assert all(sample_sheet.data.Group_By_Subject_ID == ["T0001", "T0002"])
