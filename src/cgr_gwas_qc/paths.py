from pathlib import Path
from typing import List, Sequence, Union


def make_path_list(files: Union[str, Path, Sequence[Union[str, Path]]]) -> List[Path]:
    if isinstance(files, (str, Path)):
        files = [files]

    # make sure all files are Path
    return [Path(file_) for file_ in files]
