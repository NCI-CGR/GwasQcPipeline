from pathlib import Path

import pandas as pd
import pytest

from cgr_gwas_qc import load_config
from cgr_gwas_qc.config import ConfigMgr, scan_for_yaml
from cgr_gwas_qc.testing import chdir, make_test_config
from cgr_gwas_qc.testing.data import FakeData
from cgr_gwas_qc.version import __version__


@pytest.fixture
def example_working_dir(tmp_path):
    """Returns an test working directory.

    The working directory contains:
        - Sample Sheet
        - Config File
    """
    FakeData(tmp_path).add_sample_sheet().make_config()
    return tmp_path


################################################################################
# Make sure both `config.yml` or `config.yaml` are found.
################################################################################
yaml_extensions = [
    ("config", "yml"),
    ("config", "yaml"),
]


@pytest.mark.parametrize("name,ext", yaml_extensions)
def test_scan_for_yaml(tmp_path, name, ext):
    # GIVEN: A config.yml or config.yaml
    cfg_file = tmp_path / f"{name}.{ext}"
    cfg_file.touch()

    # WHEN: we scan files in a directory
    found_config = scan_for_yaml(tmp_path, name)

    # THEN: we find the config.yml or config.yaml
    assert found_config.as_posix() == cfg_file.as_posix()


################################################################################
# Make sure if we call ConfigMgr multiple times it only creates a single
# instance.
################################################################################
def test_config_only_uses_one_instance(example_working_dir: Path, monkeypatch):
    """The ``ConfigMgr`` is normally a single-ton.

    Under normal circumstances the ``ConfigMgr`` only creates a single
    instances per python session. For testing I needed to change this
    functionality to create a new instance with each test. However, here I am
    explicitly testing this functionality so I have to undo the monkeypatch.
    """
    # GIVEN: A working directory with a config.yml and sample sheet

    # WHEN: we instantiate the ConfigMgr multiple times
    monkeypatch.undo()  # undo the ``ConfigMgr.instance()`` monkeypatch
    with chdir(example_working_dir):
        cfg = ConfigMgr.instance()
    cfg2 = ConfigMgr.instance()

    # THEN: We get the exact same object back
    assert cfg is cfg2


################################################################################
# Make sure ConfigMgr loads the yaml correctly.
################################################################################
def test_load_config_missing_sample_sheet(tmp_path):
    # GIVEN: A workding dir with a config.yml but no sample sheet
    make_test_config(tmp_path)

    # WHEN: We try to load the ConfigMgr
    # THEN: We get warnings about the missing Sample Sheet and the ``cfg.ss`` is None
    with chdir(tmp_path):
        with pytest.warns(RuntimeWarning):
            """Warn that sample_sheet.csv is not there"""
            cfg = load_config()

        with pytest.warns(RuntimeWarning):
            """Trying to access sample raises warning and gives None"""
            val = cfg.ss

        assert val is None


def test_manually_loading_config(example_working_dir: Path):
    # GIVEN: A working dir with a config.yml and a sample sheet

    # WHEN: We load the ``ConfigMgr`` using the ``instance()`` method
    with chdir(example_working_dir):
        cfg = ConfigMgr.instance()

    # THEN: We have a working ``ConfigMgr``
    assert cfg.config.pipeline_version == __version__
    assert isinstance(cfg.ss, pd.DataFrame)


def test_load_config(example_working_dir: Path):
    # GIVEN: A working dir with a config.yml and a sample sheet

    # WHEN: We load the ``ConfigMgr`` using the ``load_config()`` helper function
    with chdir(example_working_dir):
        cfg = load_config()

    # THEN: We have a working ``ConfigMgr``
    assert cfg.config.pipeline_version == __version__
    assert isinstance(cfg.ss, pd.DataFrame)
