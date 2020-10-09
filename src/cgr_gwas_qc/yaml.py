from collections import ChainMap
from pathlib import Path
from typing import List, Optional, Union

from ruamel.yaml import YAML

yaml_reader = YAML(typ="rt")


def load(files: Union[Path, List[Optional[Path]]]) -> ChainMap:
    """Load a set of Yaml configs.

    Given a list of config files loads them as a ChainMap so that values from
    the first item in the list will override values of later configs.

    Args:
        files: A config file or list of config files to be loaded. If given a
          list, ``None`` and missing files will be ignored.

    Returns:
        ChainMap: A config chian map (dict-like data structure) where values
        in the first item in the list will override latter values.
    """
    if isinstance(files, Path):
        files = [files]
    return ChainMap(*[yaml_reader.load(file.open()) for file in files if file and file.exists()])
