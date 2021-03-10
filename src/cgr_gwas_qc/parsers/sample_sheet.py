"""Parse Illumina sample sheet sections into usable pieces."""
import re
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import pandas as pd


class SampleSheet:
    """An object representation of the Illumina sample sheet.

    The Illumina sample sheet is an INI like file with different sections.
    Unlike an INI, these files are often saved as a CSV resulting in every
    row having the same number of commas. This class cleans up the sections
    by removing the extra commas, parses each section out, and converts each
    section to a more usable format.

    Attributes:
        header (dict): Key-value pairs from the header section.
        manifests (dict): Key-value pairs from the manifests section.
        data (pd.DataFrame): A dataframe representation of the data section.

    Example:
        sample_sheet = SampleSheet(sample_manifest_file_name)
    """

    section_regex = re.compile(r"(?:\[(\w+)\].*?\n([^\[]*))")

    def __init__(self, file_name: Union[str, Path]) -> None:
        """Load a sample sheet.

        Args:
            file_name: Path the an Illumina sample sheet file.
        """
        self.file_name = Path(file_name)
        self._sections = {}

        for section in self._split_sections():
            section_header, section_data = self._clean_sections(section)
            self._sections[section_header] = section_data

    @property
    def header(self) -> Dict[str, str]:
        return self._sections["header"]

    @property
    def manifests(self) -> Dict[str, str]:
        return self._sections["manifests"]

    @property
    def data(self) -> pd.DataFrame:
        return self._sections["data"]

    @data.setter
    def data(self, value: pd.DataFrame) -> None:
        self._sections["data"] = value

    def _split_sections(self) -> List[Tuple[str, str]]:
        with open(self.file_name, "r", encoding="utf-8") as fh:
            return re.findall(self.section_regex, fh.read().strip())

    def _clean_sections(
        self, section: Tuple[str, str]
    ) -> Tuple[str, Union[Dict[str, str], pd.DataFrame]]:
        header = section[0].strip().lower()

        if header == "header":
            values = self._clean_header(section[1])
        elif header == "manifests":
            values = self._clean_manifest(section[1])
        elif header == "data":
            values = self._clean_data(section[1])
        else:
            values = _strip_terminal_commas(section[1])

        return header, values

    @staticmethod
    def _clean_header(data) -> Union[Dict[str, str], str]:
        """Converts Header section into a key-value pairs."""
        stripped = _strip_terminal_commas(data)
        try:
            return _convert_to_key_value_pair(stripped)
        except Exception:
            return stripped

    @staticmethod
    def _clean_manifest(data) -> Dict[str, str]:
        """Converts Manifests section into a key-value pairs."""
        stripped = _strip_terminal_commas(data)

        res = {}
        for k, v in _convert_to_key_value_pair(stripped).items():
            if v.endswith(".bpm"):
                try:
                    # typical pattern is `snp_array`_`bpm_file`
                    res["snp_array"], res["bpm"] = v.split("_", 1)
                except ValueError:
                    # didn't follow the typical pattern to just set it as BPM
                    res["bpm"] = v
            else:
                res[k] = v

        return res

    def _clean_data(self, data) -> pd.DataFrame:
        """Converts Data section into a dataframe."""
        no_empty_rows = self._remove_empty_rows(data)
        return pd.read_csv(StringIO(no_empty_rows), low_memory=False).dropna(how="all")

    @staticmethod
    def _remove_empty_rows(data: str) -> str:
        data_no_empty_rows = re.sub("^,+$", "", data, flags=re.MULTILINE)
        data_no_extra_breaks = re.sub("\n+", "\n", data_no_empty_rows)
        return data_no_extra_breaks.lstrip()

    def add_group_by_column(self, col_name: Optional[str] = None) -> "SampleSheet":
        """Select which column in the sample sheet to use for subject grouping.

        This function adds the column `Group_By_Subject_ID` to the sample
        sheet object. The sample sheet contains multiple columns with subject
        level information. This selects which column to use for subject
        grouping. If `col_name` is provided and has a value then this will be
        used. If the sample sheet has the new `Group_By` column, then this
        column will be used if it has a value. Finally, the
        `LIMS_Individual-ID` is used as the default if neither `col_name` or
        `Group_By` are provided, or if their values are missing. Missing
        values are common for internal controls. Ideally, the `Group_By`
        column should be used in the sample sheet.

        Args:
            col_name: Typically, this is from the user config.yml
            (`workflow_params.subject_id_to-use`).

        """

        def _get_subject_id(sr: pd.Series) -> str:
            if col_name and (col_name in sr.dropna().index):
                return sr[col_name]

            if "Group_By" in sr.dropna().index:
                return sr[sr.Group_By]

            return sr["LIMS_Individual_ID"]

        self.data["Group_By_Subject_ID"] = self.data.apply(_get_subject_id, axis=1)

        return self

    def remove_Sample_IDs(self, sample_ids: Optional[Sequence[str]] = None) -> "SampleSheet":
        if sample_ids:
            self.data = self.data.query("Sample_ID not in @sample_ids")
        return self


def _strip_terminal_commas(data: str) -> str:
    """Remove all terminal commas.

    Samples sheets have a large number of commas append to the end of the
    various sections. I think this is an artifact of exporting to a CSV
    forcing each row to have the same number of columns as the main Data
    section. This strips all of these extra commas off.

    Example:
        >>> data = "key1,value1,,,,,,,,,,\nkey2,value2,,,,,,,,,,\n"
        >>> _strip_terminal_commas(data)
        "key1,value1\nkey2,value2\n
    """
    return re.sub("\n+", "\n", re.sub(r",*\n", "\n", data))


def _convert_to_key_value_pair(stripped: str) -> Dict[str, str]:
    """Converts section into key-value pairs.

    The Header and Manifests sections appear to be made up of key-value
    pairs. This converts from `key,value\n` into a dictionary.

    Args:
        stripped: A string that has all trailing commas removed.

    Example:
        >>> stripped = "key1,value1\nkey2,value2\n"
        >>> _convert_to_key_value_pair(stripped)
        {"key1": "value1", "key2": "value2"}
    """
    return {k: v for row in stripped.strip().splitlines() for k, v in [row.split(",", 1)]}
