import multiprocessing as mp
from collections import defaultdict
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from textwrap import indent
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set

import pandas as pd
import typer
from more_itertools import chunked
from pydantic.error_wrappers import ValidationError

from cgr_gwas_qc import parsers, validators, yaml
from cgr_gwas_qc.config import config_to_yaml
from cgr_gwas_qc.exceptions import GwasQcValidationError, SampleSheetNullRowError
from cgr_gwas_qc.models.config import Config, ReferenceFiles, UserFiles
from cgr_gwas_qc.reporting import CASE_CONTROL_DTYPE, SEX_DTYPE

app = typer.Typer(add_completion=False)


@dataclass(frozen=True)
class ProblemFile:
    Sample_ID: str
    reason: str
    file_type: Optional[str] = None
    filename: Optional[str] = None


@app.command()
def main(
    config_file: Path = typer.Option("config.yml", exists=True, readable=True),
    reference_check: bool = True,
    user_files_check: bool = True,
    update_config: bool = True,
    threads: int = 4,
    cluster_group_size: int = 1000,
):
    """Pre-flight checks to make sure user input files are readable and complete.

    Pre-flight checks include the Sample Sheet, reference files, Idat files,
    and GTC files. For Idat and GTC files, the user provided file name
    pattern is extended using the columns in the sample sheet.

    Pre-flight also creates a parsed version of the sample sheet
    ``cgr_sample_sheet.csv``. This file is used in the Gwas QC workflow.
    """
    config = check_config(config_file)
    ss = check_sample_sheet(
        config.sample_sheet,
        config.workflow_params.subject_id_column,
        config.workflow_params.expected_sex_column,
        config.workflow_params.case_control_column,
    )

    if reference_check:
        check_reference_files(config.reference_files)

    problem_samples = (
        check_user_files(config.user_files, ss, threads) if user_files_check else set()
    )

    if config.Sample_IDs_to_remove:
        # Remove IDs flagged in the config
        problem_samples |= {
            ProblemFile(Sample_ID, "UserExclusion") for Sample_ID in config.Sample_IDs_to_remove
        }

    # Create a parsed version of the sample sheet with some custom columns
    typer.secho("Saving Updated Sample Sheet to (cgr_sample_sheet.csv)", fg=typer.colors.GREEN)
    update_sample_sheet(
        ss,
        config.workflow_params.subject_id_column,
        config.workflow_params.expected_sex_column,
        config.workflow_params.case_control_column,
        problem_samples,
        cluster_group_size,
    ).to_csv("cgr_sample_sheet.csv", index=False)

    # Update the config file with some dynamic settings
    if update_config:
        # Update the config file with problem samples and re-calculate values
        typer.secho("Saving Updated Config to (cgr_config.yml)", fg=typer.colors.GREEN)
        update_config_file(config, ss)


def check_config(filename: Path) -> Config:
    try:
        data = yaml.load(filename)
        config = Config.parse_obj(data)
    except Exception as err:
        if isinstance(err, OSError):
            msg = err.args[1]
        elif isinstance(err, ValidationError):
            msg = str(err.args[0][0].exc).replace("\n", " ")
        else:
            msg = err.args[0]
        typer.secho(f"Config ERROR: ({filename.as_posix()})\n\t{msg}\n", fg=typer.colors.RED)
        typer.secho("Exiting... Cannot continue without a valid config file.", fg=typer.colors.RED)
        raise SystemExit

    typer.secho(f"Config OK ({filename.as_posix()})", fg=typer.colors.GREEN)
    return config


def check_reference_files(reference_files: ReferenceFiles):
    bpm_file = reference_files.illumina_manifest_file
    try:
        validators.bpm.validate(bpm_file)
        typer.secho(f"BPM OK ({bpm_file})", fg=typer.colors.GREEN)
    except Exception as err:
        msg = err.args[0] if len(err.args) == 1 else err.args[1]
        typer.secho(f"BPM ERROR: {msg} ({bpm_file})", fg=typer.colors.RED)

    vcf_file = reference_files.thousand_genome_vcf
    try:
        validators.bgzip.validate(vcf_file)
        typer.secho(f"VCF OK ({vcf_file})", fg=typer.colors.GREEN)
    except Exception as err:
        msg = err.args[0] if len(err.args) == 1 else err.args[1]
        typer.secho(f"VCF ERROR: {msg} ({vcf_file})", fg=typer.colors.RED)

    tbi_file = reference_files.thousand_genome_tbi
    try:
        validators.bgzip.validate(tbi_file)
        typer.secho(f"VCF.TBI OK ({tbi_file})", fg=typer.colors.GREEN)
    except Exception as err:
        msg = err.args[0] if len(err.args) == 1 else err.args[1]
        typer.secho(f"VCF.TBI ERROR: {msg} ({tbi_file})", fg=typer.colors.RED)


def check_sample_sheet(
    filename: Path, subject_id_column: str, expected_sex_column: str, case_control_column: str
) -> pd.DataFrame:
    try:
        if parsers.sample_sheet.is_sample_manifest(filename):
            # User provided a CGR like manifest file
            validators.sample_sheet.validate_manifest(
                filename, subject_id_column, expected_sex_column, case_control_column
            )
            df = parsers.sample_sheet.SampleManifest(filename).data
        else:
            # User provided a plain CSV file
            validators.sample_sheet.validate_sample_sheet(
                filename, subject_id_column, expected_sex_column, case_control_column
            )
            df = pd.read_csv(filename)
    except SampleSheetNullRowError:
        typer.secho(
            f"Sample Sheet WARNING: Contains Empty Rows ({filename.as_posix()})",
            fg=typer.colors.YELLOW,
        )
    except Exception as err:
        if isinstance(err, FileNotFoundError):
            msg = f"FileNotFound ({filename.as_posix()})"
        else:
            msg = err.args[1]
        typer.secho(f"Sample Sheet ERROR: {msg}", fg=typer.colors.RED)
        typer.secho("Exiting... Cannot continue without a valid sample sheet.", fg=typer.colors.RED)
        raise SystemExit

    typer.secho(f"Sample Sheet OK ({filename.as_posix()})", fg=typer.colors.GREEN)
    return df


def check_user_files(user_files: UserFiles, ss: pd.DataFrame, threads: int) -> Set[ProblemFile]:
    """Check user files (IDAT, GTC) using multiple threads.

    There can be tens of thousands of user files to process here. To speed
    things up a bit we use parallel processing. This function spins up
    ``threads`` processes and runs user files checks for each sample.

    Returns:
        Sample_IDs that had a problem with either their IDAT or GTC files.
    """
    pool = mp.Pool(threads)  # Create a pool of workers
    args = list(  # [(UserFiles, Dict[sample_sheet_column, sample_sheet_row1_value]), ...]
        product([user_files], [record._asdict() for record in ss.itertuples(index=False)])
    )
    typer.secho("Checking user files for {:,} samples.".format(len(args)))
    futures = pool.starmap(_check_user_files, args)  # Send work to pool of workers
    with typer.progressbar(futures, length=len(args)) as bar:
        problem_user_files = set()
        for results in bar:  # collect results (i.e., problem samples)
            problem_user_files |= {problem for problem in results if problem}
        _pretty_print_user_problems(problem_user_files)

    return problem_user_files


def update_config_file(config: Config, ss: pd.DataFrame):
    try:
        if config.num_snps == 0:
            config.num_snps = parsers.illumina.BeadPoolManifest(
                config.reference_files.illumina_manifest_file
            ).num_loci
    except Exception:
        typer.secho(
            "  - Problem parsing the illumina manifest file, could not update 'config.num_snp'.",
            fg=typer.colors.YELLOW,
        )

    if config.num_samples == 0:
        config.num_samples = ss.shape[0]

    config_to_yaml(config)


def update_sample_sheet(
    df: pd.DataFrame,
    subject_id_column: str,
    expected_sex_column: str,
    case_control_column: str,
    problem_samples: Iterable[ProblemFile],
    cluster_group_size: int,
) -> pd.DataFrame:
    _add_group_by_column(df, subject_id_column)
    _add_is_internal_control(df)
    _add_sample_exclusion(df, problem_samples)
    _add_user_exclusion(df, problem_samples)
    _add_missing_idats(df, problem_samples)
    _add_missing_gtc(df, problem_samples)
    _update_expected_sex(df, expected_sex_column)
    _update_case_control(df, case_control_column)
    _add_replicate_info(df)
    return _add_replicate_info(df).pipe(_add_cluster_group, cluster_group_size)


def _add_group_by_column(df: pd.DataFrame, subject_id_column: str = "Group_By"):
    """Select which column in the sample sheet to use for subject grouping.

    This function adds the column `Group_By_Subject_ID` to the sample
    sheet object. The sample sheet contains multiple columns with subject
    level information. The user defines which column to use in the config
    ``config.workflow_settings.subject_id_column``.

    Recently, Q1 2021, we started adding a column named ``Group_By`` which
    contains the column name(s) to use for subject grouping. This allows
    values from multiple columns to be used.

    Note::
        We treat ``LIMS_Individual_ID`` as the default value if other values are missing.
    """

    def _get_subject_id(sr: pd.Series) -> str:
        if "LIMS_Individual_ID" == subject_id_column:
            return sr[subject_id_column]

        default_id = sr.get("LIMS_Individual_ID", pd.NA)
        if "Group_By" == subject_id_column:
            return sr[sr[subject_id_column]] if pd.notna(sr[sr[subject_id_column]]) else default_id
        return sr[subject_id_column] if pd.notna(sr[subject_id_column]) else default_id

    df["Group_By_Subject_ID"] = df.apply(_get_subject_id, axis=1)


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


def _add_user_exclusion(df, problem_samples: Iterable[ProblemFile]):
    df["is_user_exclusion"] = False
    problem_ids = {
        problem.Sample_ID for problem in problem_samples if problem.reason == "UserExclusion"
    }
    if problem_ids:
        mask = df.Sample_ID.isin(problem_ids)
        df.loc[mask, "is_user_exclusion"] = True


def _add_missing_idats(df, problem_samples: Iterable[ProblemFile]):
    df["is_missing_idats"] = False
    problem_ids = {
        problem.Sample_ID
        for problem in problem_samples
        if problem.file_type
        and problem.file_type.startswith("idat")
        and (problem.reason == "FileNotFound")
    }
    if problem_ids:
        mask = df.Sample_ID.isin(problem_ids)
        df.loc[mask, "is_missing_idats"] = True


def _add_missing_gtc(df, problem_samples: Iterable[ProblemFile]):
    df["is_missing_gtc"] = False
    problem_ids = {
        problem.Sample_ID
        for problem in problem_samples
        if problem.file_type
        and problem.file_type.startswith("gtc")
        and (problem.reason == "FileNotFound")
    }
    if problem_ids:
        mask = df.Sample_ID.isin(problem_ids)
        df.loc[mask, "is_missing_gtc"] = True


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


def _add_sample_exclusion(df, problem_samples: Optional[Iterable[ProblemFile]]):
    df["is_sample_exclusion"] = False
    if problem_samples:
        problem_sample_ids = {problem.Sample_ID for problem in problem_samples}
        mask = df.Sample_ID.isin(problem_sample_ids)
        df.loc[mask, "is_sample_exclusion"] = True


def _add_cluster_group(df, cluster_group_size) -> pd.DataFrame:
    """Add column for grouping when running in cluster mode."""
    cluster_group = pd.Series(
        {
            v: f"cgroup{i}"
            for i, grp in enumerate(chunked(df.Sample_ID, cluster_group_size), start=1)
            for v in grp
        },
        name="cluster_group",
    ).rename_axis("Sample_ID")
    return df.merge(cluster_group, on="Sample_ID")


def _check_gtc(filename: str, sample_id: str) -> Optional[ProblemFile]:
    try:
        validators.gtc.validate(Path(filename))
    except FileNotFoundError:
        return ProblemFile(sample_id, "FileNotFound", "gtc", filename)
    except PermissionError:
        return ProblemFile(sample_id, "Permissions", "gtc", filename)
    except GwasQcValidationError as err:
        return ProblemFile(sample_id, err.args[0], "gtc", filename)
    return None  # No problems


def _check_idat(filename: str, color: str, sample_id: str) -> Optional[ProblemFile]:
    try:
        validators.idat.validate(Path(filename))
    except FileNotFoundError:
        return ProblemFile(sample_id, "FileNotFound", f"idat_{color}", filename)
    except PermissionError:
        return ProblemFile(sample_id, "Permissions", f"idat_{color}", filename)
    except GwasQcValidationError as err:
        return ProblemFile(sample_id, err.args[0], f"idat_{color}", filename)
    return None  # No problems


def _check_user_files(user_files: UserFiles, record: Mapping) -> Set[Optional[ProblemFile]]:
    """Check IDAT and GTC files for a given sample."""
    sample_id = record["Sample_ID"]
    problems = set()

    if user_files.idat_pattern:
        red_file = user_files.idat_pattern.red.format(**record)
        problems.add(_check_idat(red_file, "red", sample_id))

        green_file = user_files.idat_pattern.green.format(**record)
        problems.add(_check_idat(green_file, "green", sample_id))

    if user_files.gtc_pattern:
        gtc_file = user_files.gtc_pattern.format(**record)
        problems.add(_check_gtc(gtc_file, sample_id))

    return problems


def _filter_list(problems: Iterable[ProblemFile], file_type: str) -> Dict[str, List[str]]:
    res = defaultdict(list)
    for problem in problems:
        if problem.filename and problem.file_type == file_type:
            res[problem.reason].append(problem.filename)
    return res


def _pretty_print_paths(data: Mapping[str, Sequence[str]]) -> str:
    """For each exception output a list of files nicely."""
    output = ""
    for k, v in data.items():
        output += f"  {k}:\n"
        files = "\n".join(sorted(v))
        output += f"{indent(files, '    - ')}\n"
    return output


def _pretty_print_user_problems(problems: Iterable[ProblemFile]):
    if idat_red := _filter_list(problems, "idat_red"):
        typer.secho(
            "IDAT RED ERROR: There was a problem with these files:\n{}".format(
                _pretty_print_paths(idat_red)
            ),
            fg=typer.colors.RED,
        )
    else:
        typer.secho("IDAT RED Files OK.", fg=typer.colors.GREEN)

    if idat_green := _filter_list(problems, "idat_green"):
        typer.secho(
            "IDAT GREEN ERROR: There was a problem with these files:\n{}".format(
                _pretty_print_paths(idat_green)
            ),
            fg=typer.colors.RED,
        )
    else:
        typer.secho("IDAT GREEN Files OK.", fg=typer.colors.GREEN)

    if gtc := _filter_list(problems, "gtc"):
        typer.secho(
            "GTC ERROR: There was a problem with these files:\n{}".format(_pretty_print_paths(gtc)),
            fg=typer.colors.RED,
        )
    else:
        typer.secho("GTC Files OK.", fg=typer.colors.GREEN)


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


if __name__ == "__main__":
    app()
