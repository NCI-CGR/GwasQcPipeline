import os
import shutil
from collections import ChainMap
from pathlib import Path

import pytest

from cgr_gwas_qc.config import ConfigMgr, flatten_nested, scan_for_yaml

yaml_extensions = [
    ("config", "yml"),
    ("config", "yaml"),
    ("patterns", "yml"),
    ("patterns", "yaml"),
]


@pytest.mark.parametrize("name,ext", yaml_extensions)
def test_scan_for_yaml(tmp_path, name, ext):
    # Create a config file
    cfg_file = tmp_path / f"{name}.{ext}"
    cfg_file.touch()

    found_config = scan_for_yaml(tmp_path, name)

    assert found_config.as_posix() == cfg_file.as_posix()


flatten_items = [
    (
        {
            "lvl1": {
                "lvl1a": "lvl1a",
                "lvl1b": ["lvl1bi", "lvl1bii"],
                "lvl1c": {
                    "lvl1ca": "lvl1ca",
                    "lvl1cb": "lvl1cb",
                    "lvl1cd": ["lvl1cdi", "lvl1cdii"],
                },
            },
        },
        ["lvl1a", "lvl1bi", "lvl1bii", "lvl1ca", "lvl1cb", "lvl1cdi", "lvl1cdii"],
    ),
    (
        ChainMap({"lvl1a": "lvl1a", "lvl1b": ["lvl1bi", "lvl1bii"]}, {"lvl1c": "lvl1ca"}),
        ["lvl1a", "lvl1bi", "lvl1bii", "lvl1ca"],
    ),
]


@pytest.mark.parametrize("nested_content,expected_content", flatten_items)
def test_flatten_nested(nested_content, expected_content):
    flattened_content = flatten_nested(nested_content)
    assert sorted(flattened_content) == sorted(expected_content)


@pytest.fixture
def create_user_dir(tmpdir):
    tmp_pth = Path(tmpdir)

    # copy sample sheet
    shutil.copyfile("tests/data/example_sample_sheet.csv", tmp_pth / "samplesheet.csv")

    # create user config
    config_file = tmp_pth / "config.yml"
    config_file.write_text("project_name: test\nsample_sheet: ./samplesheet.csv")

    return tmp_pth


def test_create_config_instance(create_user_dir: Path):
    os.chdir(create_user_dir)
    cfg = ConfigMgr.instance()
    assert cfg.user_config == (create_user_dir / "config.yml").absolute()
    assert cfg.config["project_name"] == "test"

    cfg2 = ConfigMgr.instance()
    assert cfg is cfg2
