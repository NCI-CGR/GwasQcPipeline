from collections import defaultdict
from pathlib import Path
from textwrap import indent
from typing import Dict, List, Union

import typer

from cgr_gwas_qc.paths import make_path_list
from cgr_gwas_qc.validators import check_file


class IdatFileError(Exception):
    pass


def check_idat_files(files: Union[str, Path, List[Union[str, Path]]]) -> None:
    counter = defaultdict(set)
    for file_ in make_path_list(files):
        try:
            check_file(file_)
        except FileNotFoundError:
            counter["FileNotFoundError"].add(file_.as_posix())
        except PermissionError:
            counter["PermissionError"].add(file_.as_posix())

        # TODO: Add Idat file checks


def idat_summary(counter: Dict[str, set]):
    """Check the exception counter and print out issues.

    Raises:
        IdatFileError
    """
    if not counter["exception"]:
        return

    message_strings = {
        "FileNotFoundError": "Missing Idat files",
        "PermissionError": "Idat files cannot be read",
    }

    for key, files in counter.items():
        if not files:
            continue

        file_string = indent("\n".join(files), "  - ")
        typer.echo(f"{message_strings[key]} (n = {len(files):,}):\n{file_string}")

    raise IdatFileError
