from pathlib import Path

import pandas as pd
import pytest

from cgr_gwas_qc import load_config
from cgr_gwas_qc.config import ConfigMgr, scan_for_yaml
from cgr_gwas_qc.testing import chdir, make_test_config
from cgr_gwas_qc.version import __version__

################################################################################
# Make sure both `config.yml` or `config.yaml` are found.
################################################################################
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


################################################################################
# Sanity check that ConfigMgr loads the yaml correctly.
################################################################################
def test_manually_loading_config(small_bed_working_dir: Path):
    with chdir(small_bed_working_dir):
        cfg = ConfigMgr.instance()
        assert cfg.config.pipeline_version == __version__
        assert isinstance(cfg.ss, pd.DataFrame)


def test_load_config(small_bed_working_dir: Path):
    with chdir(small_bed_working_dir):
        cfg = load_config()
        assert cfg.config.pipeline_version == __version__
        assert isinstance(cfg.ss, pd.DataFrame)


def test_load_config_missing_sample_sheet(tmp_path):
    make_test_config(tmp_path)
    with chdir(tmp_path):
        with pytest.warns(RuntimeWarning):
            """Warn that sample_sheet.csv is not there"""
            cfg = load_config()

        with pytest.warns(RuntimeWarning):
            """Trying to access sample raises warning and gives None"""
            val = cfg.ss

        assert val is None


################################################################################
# Make sure if we call ConfigMgr multiple times it only creates a single
# instance.
################################################################################
def test_config_only_uses_one_instance(small_bed_working_dir: Path):
    with chdir(small_bed_working_dir):
        cfg = ConfigMgr.instance()
    cfg2 = ConfigMgr.instance()
    assert cfg is cfg2
