"""Testing parsing of the Illumina sample sheet."""
from pathlib import Path
from textwrap import dedent

import pandas as pd
import pytest

from cgr_gwas_qc.parsers.sample_sheet import (
    SampleSheet,
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
    sample_sheet = SampleSheet(file_name)

    # THEN: we apture the three section headers.
    assert sorted(sample_sheet._sections.keys()) == sorted(["header", "manifests", "data"])


################################################################################
# The parsed Header and Manifest should be dictionaries while the Data section
# should be a pandas.DataFrame.
################################################################################
def test_sample_sheet_properties_right_type(sample_sheet):
    # GIVEN: A parsed sample sheet object
    # THEN: that from each section is the correct type
    # header data is a dictionary of key-value pairs
    assert isinstance(sample_sheet.header, dict)
    # manifests data is a dictionary of key-value pairs
    assert isinstance(sample_sheet.manifests, dict)
    # data data is a dataframe
    assert isinstance(sample_sheet.data, pd.DataFrame)


################################################################################
# Sanity check that Data section gives expected results.
################################################################################
def test_sample_sheet_data(sample_sheet):
    """Test dataframe functionality."""
    # GIVEN: A parsed sample sheet object
    # THEN: the dataframe from the data section behaves as expected
    # has the same number of rows as the example sample sheet data section
    assert sample_sheet.data.shape[0] == 4
    # Allows querying by different fields and returns the right number of results
    assert sample_sheet.data.query("Identifiler_Sex == 'M'").shape[0] == 2
    assert sample_sheet.data.query("`Case/Control_Status` == 'Case'").shape[0] == 2


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
    ss = SampleSheet(sample_sheet)

    # THEN: we automatically drop that empty row
    assert ss.data.shape[0] == 2
