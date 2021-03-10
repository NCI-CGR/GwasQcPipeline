import multiprocessing as mp
from collections import defaultdict, namedtuple
from itertools import product
from pathlib import Path
from textwrap import indent
from typing import Dict, List, Mapping, Optional, Sequence, Set

import pandas as pd
import typer

from cgr_gwas_qc import load_config
from cgr_gwas_qc.config import ConfigMgr, config_to_yaml
from cgr_gwas_qc.exceptions import GwasQcValidationError
from cgr_gwas_qc.models.config import ReferenceFiles
from cgr_gwas_qc.models.config.user_files import UserFiles
from cgr_gwas_qc.parsers.illumina.IlluminaBeadArrayFiles import BeadPoolManifest
from cgr_gwas_qc.validators import bgzip, bpm, gtc, idat, sample_sheet

app = typer.Typer(add_completion=False)

ProblemFile = namedtuple("ProblemFile", "Sample_ID,reason,file_type,filename")


@app.command()
def main(
    reference_check: bool = True,
    user_files_check: bool = True,
    update_config: bool = True,
    threads: int = 4,
):
    """Pre-flight checks to make sure user input files are readable and complete.

    Pre-flight checks include the Sample Sheet, reference files, Idat files,
    and GTC files. For Idat and GTC files, the user provided file name
    pattern is extended using the columns in the sample sheet.
    """
    cfg = load_config()
    check_sample_sheet(cfg.sample_sheet_file)

    if reference_check:
        check_reference_files(cfg.config.reference_files)

    problem_samples = None
    if user_files_check:
        problem_samples = check_user_files(cfg.config.user_files, cfg.ss, threads)

    if update_config:
        # Update the config file with problem samples and re-calculate values
        update_config_file(cfg, problem_samples)


def check_sample_sheet(ss: Path):
    try:
        sample_sheet.validate(ss)
        typer.secho(f"Sample Sheet OK ({ss.as_posix()})", fg=typer.colors.GREEN)
    except sample_sheet.SampleSheetNullRowError:
        typer.secho(f"Sample Sheet Contains Empty Rows ({ss.as_posix()})", fg=typer.colors.YELLOW)


def check_reference_files(reference_files: ReferenceFiles):
    bpm_ = reference_files.illumina_manifest_file
    try:
        bpm.validate(bpm_)  # type: ignore
        typer.secho(f"BPM OK ({bpm_})", fg=typer.colors.GREEN)
    except Exception:
        typer.secho(f"BPM FAILED ({bpm_})", fg=typer.colors.RED)

    vcf_ = reference_files.thousand_genome_vcf
    try:
        bgzip.validate(vcf_)  # type: ignore
        typer.secho(f"VCF OK ({vcf_})", fg=typer.colors.GREEN)
    except Exception:
        typer.secho(f"VCF FAILED ({vcf_})", fg=typer.colors.RED)

    tbi_ = reference_files.thousand_genome_tbi
    try:
        bgzip.validate(tbi_)  # type: ignore
        typer.secho(f"VCF.TBI OK ({tbi_})", fg=typer.colors.GREEN)
    except Exception:
        typer.secho(f"VCF.TBI FAILED ({tbi_})", fg=typer.colors.RED)


def check_user_files(user_files: UserFiles, sample_sheet: pd.DataFrame, threads: int) -> Set[str]:
    # Here I use multiple processors speed up processing of user files
    pool = mp.Pool(threads)
    args = list(
        product([user_files], [record._asdict() for record in sample_sheet.itertuples(index=False)])
    )
    typer.secho("Checking user files for {:,} samples.".format(len(args)))
    futures = pool.starmap(_check_user_files, args)
    with typer.progressbar(futures, length=len(args)) as bar:
        problem_user_files = []
        for results in bar:
            problem_user_files.extend([problem for problem in results if problem])
        _pretty_print_problems(problem_user_files)

    return {problem.Sample_ID for problem in problem_user_files}


def update_config_file(cfg: ConfigMgr, problem_samples: Optional[Set[str]]):
    if problem_samples:
        # Add problem samples and update num_samples
        cfg.config.Sample_IDs_to_remove = list(problem_samples)
        cfg.config.num_samples = cfg.ss.query("Sample_ID not in @problem_samples").shape[0]

    if cfg.config.num_snps == 0:
        cfg.config.num_snps = BeadPoolManifest(
            cfg.config.reference_files.illumina_manifest_file
        ).num_loci

    config_to_yaml(cfg.config)


def _check_user_files(user_files: UserFiles, record: Mapping) -> List[Optional[ProblemFile]]:
    sample_id = record["Sample_ID"]
    problems = []

    if user_files.idat_pattern:
        problems.append(_check_idat(user_files.idat_pattern.red.format(**record), "red", sample_id))
        problems.append(
            _check_idat(user_files.idat_pattern.green.format(**record), "green", sample_id)
        )
    if user_files.gtc_pattern:
        problems.append(_check_gtc(user_files.gtc_pattern.format(**record), sample_id))

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
    return None


def _check_gtc(filename: str, sample_id: str) -> Optional[ProblemFile]:
    try:
        gtc.validate(Path(filename))
    except FileNotFoundError:
        return ProblemFile(sample_id, "FileNotFound", "gtc", filename)
    except PermissionError:
        return ProblemFile(sample_id, "Permissions", "gtc", filename)
    except GwasQcValidationError as err:
        return ProblemFile(sample_id, err.args[0], "gtc", filename)
    return None


def _pretty_print_problems(problems: Sequence[ProblemFile]):
    if idat_red := _filter_list(problems, "idat_red"):
        typer.secho(
            "There was a problem with these Idat red files:\n{}".format(
                _pretty_print_paths(idat_red)
            ),
            fg=typer.colors.RED,
        )
    else:
        typer.secho("IDAT RED Files OK.", fg=typer.colors.GREEN)

    if idat_green := _filter_list(problems, "idat_green"):
        typer.secho(
            "There was a problem with these Idat green files:\n{}".format(
                _pretty_print_paths(idat_green)
            ),
            fg=typer.colors.RED,
        )
    else:
        typer.secho("IDAT GREEN Files OK.", fg=typer.colors.GREEN)

    if gtc := _filter_list(problems, "gtc"):
        typer.secho(
            "There was a problem with these gtc files:\n{}".format(_pretty_print_paths(gtc)),
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
