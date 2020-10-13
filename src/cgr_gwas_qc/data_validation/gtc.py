from collections import defaultdict
from pathlib import Path
from textwrap import indent
from typing import Dict, List, Union

import typer

from cgr_gwas_qc.data_validation import check_file
from cgr_gwas_qc.parsers.illumina import GenotypeCalls
from cgr_gwas_qc.paths import make_path_list


class GTCFileError(Exception):
    pass


def check_gtc_files(files: Union[str, Path, List[Union[str, Path]]]) -> None:
    counter = defaultdict(set)
    for file_ in make_path_list(files):
        try:
            check_file(file_)
        except FileNotFoundError:
            counter["exception"].add(file_.as_posix())
            counter["FileNotFoundError"].add(file_.as_posix())
        except PermissionError:
            counter["exception"].add(file_.as_posix())
            counter["PermissionError"].add(file_.as_posix())

        try:
            GenotypeCalls(file_)
        except Exception as err:
            if err.args[0].startswith("GTC format error"):
                counter["exception"].add(file_.as_posix())
                counter["Not GTC File"].add(file_.as_posix())
            elif err.args[0].startswith("GTC file is incomplete"):
                counter["exception"].add(file_.as_posix())
                counter["Incomplete File"].add(file_.as_posix())
            elif err.args[0].startswith("Unsupported GTC File version"):
                counter["exception"].add(file_.as_posix())
                counter["Unsupported GTC Version"].add(file_.as_posix())

    gtc_summary(counter)


def gtc_summary(counter: Dict[str, set]):
    """Check the exception counter and print out issues.

    Raises:
        GTCFileError
    """
    if not counter["exception"]:
        return

    message_strings = {
        "FileNotFoundError": "Missing GTC Files",
        "PermissionError": "GTC files cannot be read",
        "Not GTC File": "Not a GTC file",
        "Incomplete File": "Incomplete GTC file",
        "Unsupported GTC Version": "Unsupported GTC version",
    }

    for key, files in counter.items():
        if not files:
            continue

        file_string = indent("\n".join(files), "  - ")
        typer.echo(f"{message_strings[key]} (n = {len(files):,}):\n{file_string}")

    raise GTCFileError
