from pathlib import Path
from typing import List, Union


def make_path_list(files: Union[str, Path, List[Union[str, Path]]]) -> List[Path]:
    if isinstance(files, (str, Path)):
        files = [files]

    # make sure all files are Path
    return [Path(file_) for file_ in files]


def check_file(file_name: Path):
    """A series of common checks to make sure a file is there and readable

    Raises:
        FileNoteFoundError
        PermissionError
    """

    if not file_name.exists() or not file_name.is_file():
        raise FileNotFoundError

    with file_name.open("rb") as fh:
        fh.read(1)  # raises a PermissionError if you can't read.
