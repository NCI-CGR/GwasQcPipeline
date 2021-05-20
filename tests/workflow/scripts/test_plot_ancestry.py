import pandas as pd
import pytest

from cgr_gwas_qc.workflow.scripts import plot_ancestry


def test_plot_ancestry(real_data_cache, tmp_path):
    sample_qc_csv = real_data_cache / "dev_outputs/sample_level/sample_qc.csv"
    outfile = tmp_path / "ancestry.png"
    plot_ancestry.main(sample_qc_csv, outfile)
    assert outfile.exists()


@pytest.fixture
def sample_qc_no_case(real_data_cache, tmp_path):
    sample_qc_csv = real_data_cache / "dev_outputs/sample_level/sample_qc.csv"
    outfile = tmp_path / "sample_qc.csv"
    pd.read_csv(sample_qc_csv).query("case_control != 'Case'").to_csv(outfile, index=False)
    return outfile


def test_plot_ancestry_no_case(sample_qc_no_case, tmp_path):
    outfile = tmp_path / "ancestry.png"
    plot_ancestry.main(sample_qc_no_case, outfile)
    assert outfile.exists()


@pytest.fixture
def sample_qc_no_control(real_data_cache, tmp_path):
    sample_qc_csv = real_data_cache / "dev_outputs/sample_level/sample_qc.csv"
    outfile = tmp_path / "sample_qc.csv"
    pd.read_csv(sample_qc_csv).query("case_control != 'Control'").to_csv(outfile, index=False)
    return outfile


def test_plot_ancestry_no_controls(sample_qc_no_control, tmp_path):
    outfile = tmp_path / "ancestry.png"
    plot_ancestry.main(sample_qc_no_control, outfile)
    assert outfile.exists()
