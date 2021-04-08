from pathlib import Path
from textwrap import dedent

import pandas as pd
import pytest

from cgr_gwas_qc import load_config
from cgr_gwas_qc.config import ConfigMgr, scan_for_yaml
from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.testing.data import FakeData
from cgr_gwas_qc.version import __version__


@pytest.fixture
def example_working_dir(tmp_path):
    """Returns an test working directory.

    The working directory contains:
        - Sample Sheet
        - Config File
    """
    FakeData(tmp_path).copy_sample_sheet().make_config()
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


################################################################################
# Assign grouping column (see Issue #32)
################################################################################
def test_group_by_column_default(tmp_path):
    """Populate the Group_By_Subject_ID column with default

    The default value should be the same as LIMS_Individual_ID.
    """
    # GIVEN: Fake sample sheet and config
    FakeData(tmp_path).copy_sample_sheet().make_config()

    # WHEN: I load the config and sample sheet
    with chdir(tmp_path):
        cfg = load_config()

    # THEN: The `Group_By_Subject_ID` column is populated by the
    # `LIMS_Individual_ID` by default.
    assert all(cfg.ss.Group_By_Subject_ID == cfg.ss.LIMS_Individual_ID)


def test_group_by_column_config_option(tmp_path):
    """Populate Group_By_Subject_ID column with config option.

    The value should be the same as the column specified in
    `workflow_params.subject_id_column`.
    """
    # GIVEN: Fake sample sheet and a config where I set the `subject_id_column`
    FakeData(tmp_path).copy_sample_sheet().make_config(
        workflow_params=dict(subject_id_column="Sample_ID")
    )

    # WHEN: I load the config and sample sheet
    with chdir(tmp_path):
        cfg = load_config()

    # THEN: The `Group_By_Subject_ID` column is populated by the
    # which ever column is specified here.
    assert all(cfg.ss.Group_By_Subject_ID == cfg.ss.Sample_ID)


def test_group_by_column(tmp_path):
    """Populate Group_By_Subject_ID column with config option.

    The value should be the same as the column specified in
    `workflow_params.subject_id_column`.
    """
    # GIVEN: Fake sample sheet and config
    FakeData(tmp_path).copy_sample_sheet().make_config()
    # and the sample sheet has the `Group_By` column set
    (tmp_path / "sample_sheet.csv").write_text(
        dedent(
            """\
            [Header],,
            test,data,
            [Manifests],,
            test,data,
            [Data]
            Sample_ID,LIMS_Individual_ID,Group_By
            T0001,L0001,Sample_ID
            T0002,L0002,LIMS_Individual_ID
            """
        )
    )

    # WHEN: I load the config and sample sheet
    with chdir(tmp_path):
        cfg = load_config()

    # THEN: The `Group_By_Subject_ID` column is populated by the
    # which ever column is specified for each sample separately.
    assert all(cfg.ss.Group_By_Subject_ID == ["T0001", "L0002"])
