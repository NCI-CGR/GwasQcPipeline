from collections import defaultdict
from pathlib import Path
from textwrap import indent
from typing import Dict, List

import typer
from snakemake.rules import expand

from cgr_gwas_qc import load_config
from cgr_gwas_qc.exceptions import GwasQcValidationError
from cgr_gwas_qc.models.config import ReferenceFiles
from cgr_gwas_qc.validators import bgzip, bpm, gtc

app = typer.Typer(add_completion=False)


@app.command()
def main(reference_check: bool = True, gtc_check: bool = True, idat_check: bool = False):
    """Pre-flight checks to make sure user input files are readable and complete.

    Pre-flight checks include the Sample Sheet, reference files, Idat files,
    and GTC files. For Idat and GTC files, the user provided file name
    pattern is extended using the columns in the sample sheet.
    """
    cfg = load_config()  # This already validates the sample sheet
    typer.secho(f"Sample Sheet OK ({cfg.config.sample_sheet})", fg=typer.colors.GREEN)
    if reference_check:
        check_reference_files(cfg.config.reference_paths)

    if idat_check:
        typer.secho("Idat pre-flight is not impelmented yet.", fg=typer.colors.YELLOW)
        # check_idat_files(expand(cfg.config.user_data_patterns.idat.red, **cfg.ss.to_dict("list")), "Red")
        # check_idat_files(expand(cfg.config.user_data_patterns.idat.grn, **cfg.ss.to_dict("list")), "Green")

    if gtc_check:
        check_gtc_files(expand(cfg.config.user_data_patterns.gtc, zip, **cfg.ss.to_dict("list")))


def check_reference_files(reference_paths: ReferenceFiles):
    bpm_ = reference_paths.illumina_array_manifest
    try:
        bpm.validate(bpm_)
        typer.secho(f"BPM OK ({bpm_})", fg=typer.colors.GREEN)
    except Exception:
        typer.secho(f"BPM FAILED ({bpm_})", fg=typer.colors.RED)

    vcf_ = reference_paths.thousand_genome_vcf
    try:
        bgzip.validate(vcf_)
        typer.secho(f"VCF OK ({vcf_})", fg=typer.colors.GREEN)
    except Exception:
        typer.secho(f"VCF FAILED ({vcf_})", fg=typer.colors.RED)

    tbi_ = reference_paths.thousand_genome_tbi
    try:
        bgzip.validate(tbi_)
        typer.secho(f"VCF.TBI OK ({tbi_})", fg=typer.colors.GREEN)
    except Exception:
        typer.secho(f"VCF.TBI FAILED ({tbi_})", fg=typer.colors.RED)


def check_idat_files(files: List[str], color: str):
    problems = defaultdict(list)
    for idat in files:
        try:
            raise NotImplementedError
        except FileNotFoundError:
            problems["FileNotFound"].append(idat)
        except PermissionError:
            problems["FileNotReadable"].append(idat)
        except GwasQcValidationError as err:
            problems[err.args[0]].append(idat)

    if problems:
        typer.echo(
            "There was a problem with these Idat {} files:\n{}".format(
                color, pretty_print_paths(problems)
            )
        )
    else:
        typer.echo(f"Idat {color} Files OK.")


def check_gtc_files(files: List[str]):
    problems = defaultdict(list)
    typer.echo("Processing GTC Files")
    with typer.progressbar(files) as progress:
        for gtc_ in progress:
            try:
                gtc.validate(Path(gtc_))
            except FileNotFoundError:
                problems["FileNotFound"].append(gtc_)
            except PermissionError:
                problems["FileNotReadable"].append(gtc_)
            except GwasQcValidationError as err:
                problems[err.args[0]].append(gtc_)

    if problems:
        typer.secho(
            "There was a problem with these GTC files:\n{}".format(pretty_print_paths(problems)),
            fg=typer.colors.RED,
        )
    else:
        typer.secho(f"{len(files):,} GTC Files OK.", fg=typer.colors.GREEN)


def pretty_print_paths(data: Dict[str, List[str]]) -> str:
    """For each exception output a list of files nicely."""
    output = ""
    for k, v in data.items():
        output += f"  {k}:\n"
        files = "\n".join(sorted(v))
        output += f"{indent(files, '    - ')}\n"
    return output


if __name__ == "__main__":
    app()
