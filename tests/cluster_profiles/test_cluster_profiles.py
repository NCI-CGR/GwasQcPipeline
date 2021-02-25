import shutil
from pathlib import Path

import pytest

from cgr_gwas_qc.cluster_profiles import (
    _get_rule_names,
    load_cluster_config,
    update_group_properties,
    update_properties,
)


@pytest.mark.parametrize(
    "config_string,config",
    [
        pytest.param(
            "__default__:\n  mem_gb: 4", {"__default__": {"mem_gb": 4}}, id="basic_config"
        ),
        pytest.param(
            "__resources__:\n  mem_gb: 4", {"__default__": {}}, id="no __default__ section"
        ),
    ],
)
def test_load_cluster_config(config_string, config, tmp_path):
    # GIVEN: A cluster config file
    config_file = tmp_path / "config.yml"
    config_file.write_text(config_string)

    # WHEN: load config as dictionary
    res = load_cluster_config(config_file)

    # THEN: match expected dictionary
    assert config["__default__"] == res["__default__"]


@pytest.mark.parametrize(
    "cluster_config,job_properties,expected",
    [
        pytest.param(
            {},
            {"rule": "test_rule", "jobid": 1, "threads": 1},
            {"rulename": "test_rule", "job_id": 1, "threads": 1},
            id="basic properties",
        ),
        pytest.param(
            {"test_rule": {"threads": 2}},
            {"rule": "test_rule", "jobid": 1, "threads": 1},
            {"rulename": "test_rule", "job_id": 1, "threads": 2},
            id="cluster profile adds threads",
        ),
        pytest.param(
            {},
            {"rule": "test_rule", "jobid": 1, "threads": 1, "resources": {"mem_gb": 4}},
            {"rulename": "test_rule", "job_id": 1, "threads": 1, "mem_gb": 4},
            id="mem resources added",
        ),
        pytest.param(
            {},
            {
                "rule": "test_rule",
                "jobid": 1,
                "threads": 1,
                "resources": {"mem_gb": 4},
                "cluster": {"mem_gb": 1},
            },
            {"rulename": "test_rule", "job_id": 1, "threads": 1, "mem_gb": 4},
            id="resources override --cluster-config",
        ),
    ],
)
def test_update_properties(cluster_config, job_properties, expected):
    options = {}
    update_properties(options, cluster_config, job_properties)
    assert expected == options


@pytest.fixture
def basic_group(tmp_path) -> Path:
    pth: Path = tmp_path / "job_script.sh"
    shutil.copyfile("tests/data/job_scripts/basic_group.sh", pth)
    return pth


@pytest.fixture
def multisample_group(tmp_path) -> Path:
    pth: Path = tmp_path / "job_script.sh"
    shutil.copyfile("tests/data/job_scripts/single_rule_10_samples_j1.sh", pth)
    return pth


def test_get_rule_names_basic_group(basic_group):
    rule_names = set(_get_rule_names(basic_group))
    assert set(["a", "b"]) == rule_names


def test_get_rule_names_multisample_group(multisample_group):
    rule_names = set(_get_rule_names(multisample_group))
    assert set(["a"]) == rule_names


@pytest.mark.parametrize(
    "cluster_config,job_properties,expected",
    [
        pytest.param(
            {},
            {"type": "group", "jobid": 1, "threads": 1},
            {"rulename": "GROUP.a.b", "job_id": 1, "threads": 1},
            id="basic group",
        ),
        pytest.param(
            {},
            {"type": "group", "jobid": 1, "threads": 1, "resources": {"mem_gb": 2}},
            {"rulename": "GROUP.a.b", "job_id": 1, "threads": 1, "mem_gb": 2},
            id="basic group with resources",
        ),
    ],
)
def test_update_group_properties_basic_group(cluster_config, job_properties, expected, basic_group):
    options = {}
    update_group_properties(options, cluster_config, job_properties, basic_group)
    assert expected == options
