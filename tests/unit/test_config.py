from pathlib import Path

import pandas as pd
import pytest

from cgr_gwas_qc.config import ConfigMgr, scan_for_yaml
from cgr_gwas_qc.testing import chdir

yaml_extensions = [
    ("config", "yml"),
    ("config", "yaml"),
]


@pytest.mark.parametrize("name,ext", yaml_extensions)
def test_scan_for_yaml(tmp_path, name, ext):
    # Create a config file
    cfg_file = tmp_path / f"{name}.{ext}"
    cfg_file.touch()

    found_config = scan_for_yaml(tmp_path, name)

    assert found_config.as_posix() == cfg_file.as_posix()


def test_config_loads(working_dir: Path):
    with chdir(working_dir):
        cfg = ConfigMgr.instance()
        assert cfg.user_config == (working_dir / "config.yml").absolute()
        assert cfg.config.project_name == "Test Project"
        assert isinstance(cfg.ss, pd.DataFrame)


def test_config_only_uses_one_instance(working_dir: Path):
    with chdir(working_dir):
        cfg = ConfigMgr.instance()
    cfg2 = ConfigMgr.instance()
    assert cfg is cfg2
