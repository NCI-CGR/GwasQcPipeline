import multiprocessing as mp
from collections import defaultdict
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from textwrap import indent
from typing import Dict, List, Mapping, Optional, Sequence, Set

import pandas as pd
import typer
from pydantic.error_wrappers import ValidationError

from cgr_gwas_qc import yaml
from cgr_gwas_qc.config import config_to_yaml
from cgr_gwas_qc.exceptions import GwasQcValidationError, SampleSheetNullRowError
from cgr_gwas_qc.models.config import Config, ReferenceFiles, UserFiles
from cgr_gwas_qc.parsers.illumina.IlluminaBeadArrayFiles import BeadPoolManifest
from cgr_gwas_qc.parsers.sample_sheet import SampleManifest, is_sample_manifest, update_sample_sheet
from cgr_gwas_qc.validators import bgzip, bpm, gtc, idat, sample_sheet

app = typer.Typer(add_completion=False)


@dataclass
class ProblemFile:
    Sample_ID: str
    reason: str
    file_type: str
    filename: str


@app.command()
def main(
    config_file: Path = typer.Option("config.yml", exists=True, readable=True),
    reference_check: bool = True,
    user_files_check: bool = True,
    update_config: bool = True,
    threads: int = 4,
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
        problem_samples |= set(config.Sample_IDs_to_remove)

    # Create a parsed version of the sample sheet with some custom columns
    typer.secho("Saving Updated Sample Sheet to (cgr_sample_sheet.csv)", fg=typer.colors.GREEN)
    update_sample_sheet(
        ss,
        config.workflow_params.subject_id_column,
        config.workflow_params.expected_sex_column,
        config.workflow_params.case_control_column,
        problem_samples,
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


def check_sample_sheet(
    filename: Path, subject_id_column: str, expected_sex_column: str, case_control_column: str
) -> pd.DataFrame:
    try:
        if is_sample_manifest(filename):
            # User provided a CGR like manifest file
            sample_sheet.validate_manifest(
                filename, subject_id_column, expected_sex_column, case_control_column
            )
            df = SampleManifest(filename).data
        else:
            # User provided a plain CSV file
            sample_sheet.validate_sample_sheet(
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


def check_reference_files(reference_files: ReferenceFiles):
    bpm_file = reference_files.illumina_manifest_file
    try:
        bpm.validate(bpm_file)
        typer.secho(f"BPM OK ({bpm_file})", fg=typer.colors.GREEN)
    except Exception as err:
        msg = err.args[0] if len(err.args) == 1 else err.args[1]
        typer.secho(f"BPM ERROR: {msg} ({bpm_file})", fg=typer.colors.RED)

    vcf_file = reference_files.thousand_genome_vcf
    try:
        bgzip.validate(vcf_file)
        typer.secho(f"VCF OK ({vcf_file})", fg=typer.colors.GREEN)
    except Exception as err:
        msg = err.args[0] if len(err.args) == 1 else err.args[1]
        typer.secho(f"VCF ERROR: {msg} ({vcf_file})", fg=typer.colors.RED)

    tbi_file = reference_files.thousand_genome_tbi
    try:
        bgzip.validate(tbi_file)
        typer.secho(f"VCF.TBI OK ({tbi_file})", fg=typer.colors.GREEN)
    except Exception as err:
        msg = err.args[0] if len(err.args) == 1 else err.args[1]
        typer.secho(f"VCF.TBI ERROR: {msg} ({tbi_file})", fg=typer.colors.RED)


def check_user_files(user_files: UserFiles, ss: pd.DataFrame, threads: int) -> Set[str]:
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
        problem_user_files = []
        for results in bar:  # collect results (i.e., problem samples)
            problem_user_files.extend([problem for problem in results if problem])
        _pretty_print_user_problems(problem_user_files)

    return {problem.Sample_ID for problem in problem_user_files}


def update_config_file(config: Config, ss: pd.DataFrame):
    try:
        if config.num_snps == 0:
            config.num_snps = BeadPoolManifest(
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


def _check_user_files(user_files: UserFiles, record: Mapping) -> List[Optional[ProblemFile]]:
    """Check IDAT and GTC files for a given sample."""
    sample_id = record["Sample_ID"]
    problems = []

    if user_files.idat_pattern:
        red_file = user_files.idat_pattern.red.format(**record)
        problems.append(_check_idat(red_file, "red", sample_id))

        green_file = user_files.idat_pattern.green.format(**record)
        problems.append(_check_idat(green_file, "green", sample_id))

    if user_files.gtc_pattern:
        gtc_file = user_files.gtc_pattern.format(**record)
        problems.append(_check_gtc(gtc_file, sample_id))

    return problems


def _check_idat(filename: str, color: str, sample_id: str) -> Optional[ProblemFile]:
    try:
        idat.validate(Path(filename))
    except FileNotFoundError:
        return ProblemFile(sample_id, "FileNotFound", f"idat_{color}", filename)
    except PermissionError:
        return ProblemFile(sample_id, "Permissions", f"idat_{color}", filename)
    except GwasQcValidationError as err:
        return ProblemFile(sample_id, err.args[0], f"idat_{color}", filename)
    return None  # No problems


def _check_gtc(filename: str, sample_id: str) -> Optional[ProblemFile]:
    try:
        gtc.validate(Path(filename))
    except FileNotFoundError:
        return ProblemFile(sample_id, "FileNotFound", "gtc", filename)
    except PermissionError:
        return ProblemFile(sample_id, "Permissions", "gtc", filename)
    except GwasQcValidationError as err:
        return ProblemFile(sample_id, err.args[0], "gtc", filename)
    return None  # No problems


def _pretty_print_user_problems(problems: Sequence[ProblemFile]):
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


def _filter_list(problems: Sequence[ProblemFile], file_type: str) -> Dict[str, List[str]]:
    res = defaultdict(list)
    for problem in problems:
        if problem.file_type == file_type:
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


if __name__ == "__main__":
    app()
