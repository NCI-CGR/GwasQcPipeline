from collections import defaultdict
from logging import getLogger
from pathlib import Path
from textwrap import indent
from typing import Dict, List

import typer
from snakemake.rules import expand

from cgr_gwas_qc import load_config
from cgr_gwas_qc.validators import GwasQcValidationError
from cgr_gwas_qc.validators.gtc import validate as gtc_validate

app = typer.Typer(add_completion=False)
logger = getLogger(__name__)


@app.command()
def main(gtc_check: bool = True, idat_check: bool = False):
    """Pre-flight checks to make sure user input files are readable and complete.

    Pre-flight checks include the Sample Sheet, reference files, Idat files,
    and GTC files. For Idat and GTC files, the user provided file name
    pattern is extended using the columns in the sample sheet.
    """
    cfg = load_config()  # This already validates all Reference files in the config
    typer.echo("Sample Sheet and Reference Files OK.")

    if idat_check:
        logger.warn("Idat pre-flight is not impelmented yet.")
        # check_idat_files(expand(cfg.config.user_data_patterns.idat.red, **cfg.ss.to_dict("list")), "Red")
        # check_idat_files(expand(cfg.config.user_data_patterns.idat.grn, **cfg.ss.to_dict("list")), "Green")

    if gtc_check:
        check_gtc_files(expand(cfg.config.user_data_patterns.gtc, **cfg.ss.to_dict("list")))


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
    for gtc in files:
        try:
            gtc_validate(Path(gtc))
        except FileNotFoundError:
            problems["FileNotFound"].append(gtc)
        except PermissionError:
            problems["FileNotReadable"].append(gtc)
        except GwasQcValidationError as err:
            problems[err.args[0]].append(gtc)

    if problems:
        typer.echo(
            "There was a problem with these GTC files:\n{}".format(pretty_print_paths(problems))
        )
    else:
        typer.echo("GTC Files OK.")


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
