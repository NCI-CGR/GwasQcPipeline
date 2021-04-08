"""Parse Illumina sample sheet sections into usable pieces."""
import re
from io import StringIO
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Union

import pandas as pd

from cgr_gwas_qc.reporting.constants import CASE_CONTROL_DTYPE, SEX_DTYPE
from cgr_gwas_qc.typing import PathLike

DTYPES = {
    "Sample_ID": "string",
    "Group_By_Subject_ID": "string",
    "num_samples_per_subject": "UInt8",
    "replicate_ids": "string",
    "expected_sex": SEX_DTYPE,
    "case_control": CASE_CONTROL_DTYPE,
    "is_internal_control": "boolean",
    "is_sample_exclusion": "boolean",
}


def read(
    filename: PathLike, remove_exclusions: bool = True, all_user_column: bool = False
) -> pd.DataFrame:
    df = pd.read_csv(filename, dtype=DTYPES)

    if not all_user_column:
        df = df.reindex(DTYPES.keys(), axis=1)

    if remove_exclusions:
        df = df.query("not is_sample_exclusion")

    return df


class SampleManifest:
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

    def __init__(self, filename: PathLike) -> None:
        """Load a sample sheet.

        Args:
            file_name: Path the an Illumina sample sheet file.
        """
        self.file_name = Path(filename)
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


def is_sample_manifest(filename: Path) -> bool:
    text = filename.read_text()
    return all(["[Header]" in text, "[Manifests]" in text, "[Data]" in text])


def update_sample_sheet(
    df: pd.DataFrame,
    subject_id_column: str,
    expected_sex_column: str,
    case_control_column: str,
    problem_sample_ids: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    _add_group_by_column(df, subject_id_column)
    _add_is_internal_control(df)
    _add_sample_exclusion(df, problem_sample_ids)
    _update_expected_sex(df, expected_sex_column)
    _update_case_control(df, case_control_column)
    return _add_replicate_info(df)


def _add_group_by_column(df: pd.DataFrame, subject_id_column: str = "Group_By"):
    """Select which column in the sample sheet to use for subject grouping.

    This function adds the column `Group_By_Subject_ID` to the sample
    sheet object. The sample sheet contains multiple columns with subject
    level information. The user defines which column to use in the config
    ``config.workflow_settings.subject_id_column``.

    Recently, Q1 2021, we started adding a column named ``Group_By`` which
    contains the column name(s) to use for subject grouping. This allows
    values from multiple columns to be used.
    """

    def _get_subject_id(sr: pd.Series) -> str:
        if "Group_By" == subject_id_column:
            return sr[sr[subject_id_column]]

        return sr[subject_id_column]

    df["Group_By_Subject_ID"] = df.apply(_get_subject_id, axis=1)


def _add_replicate_info(df: pd.DataFrame) -> pd.DataFrame:
    """Adds information about the number of samples per subject.

    This function adds two columns:
        - ``num_samples_per_subject`` count of the number of samples per subject
        - ``replicate_ids`` Concatenated Sample_IDs for samples from the same subject
    """

    def _replicate_info(x):
        num_reps = x.shape[0]
        x["num_samples_per_subject"] = num_reps
        x["replicate_ids"] = pd.NA

        if num_reps > 1:
            x["replicate_ids"] = x.Sample_ID.sort_values().str.cat(sep="|")

        return x

    return df.groupby("Group_By_Subject_ID", dropna=False).apply(_replicate_info)


def _add_sample_exclusion(df, problem_sample_ids: Optional[Iterable[str]]):
    df["is_sample_exclusion"] = False
    if problem_sample_ids:
        mask = df.Sample_ID.isin(problem_sample_ids)
        df.loc[mask, "is_sample_exclusion"] = True


def _add_is_internal_control(df: pd.DataFrame):
    """Add a flag if a sample is really an internal control.

    This function adds the column ``is_internal_control`` if it does not
    already exist. This column is generated based on assumptions about CGRs
    LIMS sheet. Third party users should create this column separately if
    they have their own internal controls.
    """
    if "is_internal_control" in df.columns:
        # Column already exists, this would allow people to add their own if they want
        return

    if "Sample_Group" in df.columns:
        df["is_internal_control"] = (df.Sample_Group == "sVALD-001").astype("boolean")
        return

    df["is_internal_control"] = False


def _update_expected_sex(df: pd.DataFrame, expected_sex_column: str = "Expected_Sex"):
    """Normalize sex calls.

    This function creates a new column ``expected_sex`` based on the column
    provided by the user ``config.workflow_settings.expected_sex_column``.
    Here we normalize sex values and update sex for internal controls.
    """
    sex_mapper = {"m": "M", "male": "M", "f": "F", "female": "F"}
    df["expected_sex"] = (
        df[expected_sex_column]
        .str.lower()
        .map(lambda sex: sex_mapper.get(sex, "U"))
        .astype(SEX_DTYPE)
    )

    if "Identifiler_Sex" in df.columns:
        # For internal controls use the `Indentifiler_Sex` column as `expected_sex`
        df.loc[df.is_internal_control, "expected_sex"] = df.loc[
            df.is_internal_control, "Identifiler_Sex"
        ]


def _update_case_control(df: pd.DataFrame, case_control_column: str = "Case/Control_Status"):
    """Normalize Case/Control annotations.

    This function creates a new column ``case_control`` based on the column
    provided by the user ``config.workflow_settings.case_control_column``.
    Here we normalize labels values and update labels to ``QC`` for internal
    controls.
    """
    case_control_mapper = {cat.lower(): cat for cat in CASE_CONTROL_DTYPE.categories}
    df["case_control"] = (
        df[case_control_column]
        .str.lower()
        .map(lambda status: case_control_mapper.get(status, "Unknown"))
        .astype(CASE_CONTROL_DTYPE)
    )

    # For internal controls set case_control to qc
    df.loc[df.is_internal_control, "case_control"] = case_control_mapper["qc"]


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
