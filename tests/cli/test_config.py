from cgr_gwas_qc.cli.config import CGEMS_QC_DIR, _get_cgems_production_run_dir


def test_get_cgems_production_run_dir(monkeypatch):
    monkeypatch.setattr("cgr_gwas_qc.cli.config.TODAY", "01012021")
    test_path = _get_cgems_production_run_dir("TestProject")
    assert CGEMS_QC_DIR / "TestProject/builds/QC_v1_01012021" == test_path


def test_get_cgems_production_run_dir_increment_if_previous_runs(monkeypatch):
    def glob(*args, **kwargs):
        """Return 2 previous runs"""
        return ["QC_v1", "QC_v2"]

    monkeypatch.setattr("cgr_gwas_qc.cli.config.Path.glob", glob)
    monkeypatch.setattr("cgr_gwas_qc.cli.config.TODAY", "01012021")
    test_path = _get_cgems_production_run_dir("TestProject")
    assert CGEMS_QC_DIR / "TestProject/builds/QC_v3_01012021" == test_path
