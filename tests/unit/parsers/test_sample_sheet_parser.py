"""Testing parsing of the Illumina sample sheet."""
from pathlib import Path

import pandas as pd
import pytest

from cgr_gwas_qc.parsers.sample_sheet import (
    SampleSheet,
    _convert_to_key_value_pair,
    _strip_terminal_commas,
)

strings_with_extra_commas = [
    ("key1,value1,,,,,,,,,,\nkey2,value2,,,,,,,,,,\n", "key1,value1\nkey2,value2\n"),
    ("key1,value1,,,,,,,,,,\n,,,,,,,,,,,\n,,,,,,,,,,,\n", "key1,value1\n"),
]


@pytest.mark.parametrize("data,expected", strings_with_extra_commas)
def test_strip_terminal_commas(data, expected):
    assert _strip_terminal_commas(data) == expected


strings_without_extra_commas = [
    ("key1,value1\nkey2,value2\n", {"key1": "value1", "key2": "value2"}),
    ("key1,value1\n", {"key1": "value1"}),
]


@pytest.mark.parametrize("data,expected", strings_without_extra_commas)
def test_convert_to_key_value_pair(data, expected):
    assert _convert_to_key_value_pair(data) == expected


example_sample_sheet = "tests/data/example_sample_sheet.csv"


@pytest.mark.parametrize("file_name", [example_sample_sheet, Path(example_sample_sheet)])
def test_load_str_or_path(file_name):
    """Make sure SampleSheet works with strings or paths."""
    _sample_sheet = SampleSheet(file_name)
    assert sorted(_sample_sheet._sections.keys()) == sorted(["header", "manifests", "data"])


@pytest.fixture(scope="module")
def sample_sheet(sample_sheet_file):
    return SampleSheet(sample_sheet_file)


def test_sample_sheet_properties_right_type(sample_sheet):
    assert isinstance(sample_sheet.header, dict)
    assert isinstance(sample_sheet.manifests, dict)
    assert isinstance(sample_sheet.data, pd.DataFrame)


def test_sample_sheet_data(sample_sheet):
    """Test dataframe functionality."""
    assert sample_sheet.data.shape[0] == 4
    assert sample_sheet.data.query("Identifiler_Sex == 'M'").shape[0] == 2
    assert sample_sheet.data.query("`Case/Control_Status` == 'Case'").shape[0] == 2
