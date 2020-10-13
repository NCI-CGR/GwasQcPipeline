from pathlib import Path


def check_file(file_name: Path) -> None:
    """A series of common checks to make sure a file is there and readable

    Raises:
        FileNoteFoundError
        PermissionError
    """

    if not file_name.exists() or not file_name.is_file():
        raise FileNotFoundError

    with file_name.open("rb") as fh:
        fh.read(1)  # raises a PermissionError if you can't read.
