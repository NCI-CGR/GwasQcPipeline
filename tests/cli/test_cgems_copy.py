from pathlib import Path

import pytest

from cgr_gwas_qc.cli import cgems_copy


@pytest.fixture(scope="module")
def gsa_test_path(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("gsa_lab_test")
    (tmp_path / "TestProject/builds/QC_v1_01012021").mkdir(parents=True, exist_ok=True)
    (tmp_path / "TestProject/builds/QC_v2_02022021").mkdir(parents=True, exist_ok=True)
    (tmp_path / "TestProject/builds/QC_v2_02022021/config.yml").touch()
    (tmp_path / "TestProject/builds/QC_v2_02022021/cgr_sample_sheet.csv").touch()
    return tmp_path


@pytest.mark.parametrize(
    "path",
    [
        pytest.param("TestProject", id="project"),
        pytest.param("TestProject/builds", id="build"),
        pytest.param("TestProject/builds/QC_v1_01012021", id="run"),
    ],
)
def test_project_name(gsa_test_path, path):
    assert "TestProject" == cgems_copy._get_project_name(gsa_test_path / path)


def test_get_previous_run_dir(gsa_test_path, monkeypatch):
    monkeypatch.setattr("cgr_gwas_qc.cli.cgems_copy.CGEMS_QC_DIR", gsa_test_path)
    assert Path(
        gsa_test_path / "TestProject/builds/QC_v2_02022021"
    ) == cgems_copy._get_previous_cgems_production_run_dir("TestProject")


def test_clone_previous_run(gsa_test_path, monkeypatch):
    monkeypatch.setattr("cgr_gwas_qc.cli.cgems_copy.CGEMS_QC_DIR", gsa_test_path)
    monkeypatch.setattr("cgr_gwas_qc.cli.config.CGEMS_QC_DIR", gsa_test_path)
    monkeypatch.setattr("cgr_gwas_qc.cli.config.TODAY", "03032021")

    project_name = "TestProject"
    previous_run_dir = cgems_copy._get_previous_cgems_production_run_dir(project_name)
    run_dir = cgems_copy._get_cgems_production_run_dir(project_name)
    cgems_copy._create_cgems_production_run_dir(run_dir, no_prompt=True)
    cgems_copy._clone_previous_run(previous_run_dir, run_dir)

    assert gsa_test_path / "TestProject/builds/QC_v3_03032021" == run_dir
    assert (run_dir / "config.yml").exists()
    assert (run_dir / "cgr_sample_sheet.csv").exists()
