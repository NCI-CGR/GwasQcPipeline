from pathlib import Path
from typing import List, Union


def make_path_list(files: Union[str, Path, List[Union[str, Path]]]) -> List[Path]:
    if isinstance(files, (str, Path)):
        files = [files]

    # make sure all files are Path
    return [Path(file_) for file_ in files]
