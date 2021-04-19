from cgr_gwas_qc.workflow.scripts import plot_ancestry


def test_plot_ancestry(sample_qc_csv, tmp_path):
    outfile = tmp_path / "ancestry.png"
    plot_ancestry.main(sample_qc_csv, outfile)
    assert outfile.exists()


def test_plot_ancestry_no_case(sample_qc_df, tmp_path):
    sample_qc_csv = tmp_path / "qc.csv"
    outfile = tmp_path / "ancestry.png"
    sample_qc_df.query("case_control != 'Case'").to_csv(sample_qc_csv, index=False)
    plot_ancestry.main(sample_qc_csv, outfile)
    assert outfile.exists()


def test_plot_ancestry_no_controls(sample_qc_df, tmp_path):
    sample_qc_csv = tmp_path / "qc.csv"
    outfile = tmp_path / "ancestry.png"
    sample_qc_df.query("case_control != 'Control'").to_csv(sample_qc_csv, index=False)
    plot_ancestry.main(sample_qc_csv, outfile)
    assert outfile.exists()
