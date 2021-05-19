import pytest

from cgr_gwas_qc.workflow.scripts import plot_autosomal_heterozygosity


@pytest.mark.real_data
def test_plot(real_data_cache, software_params, tmp_path):
    subject_qc = real_data_cache / "dev_outputs/subject_level/subjects_for_analysis.csv"
    het = real_data_cache / "dev_outputs/subject_level/European/subjects.het"
    outfile = tmp_path / "test.png"
    plot_autosomal_heterozygosity.main(
        subject_qc, het, "European", software_params.autosomal_het_threshold, outfile
    )

    assert outfile.exists()
