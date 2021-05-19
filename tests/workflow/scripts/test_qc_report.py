import pytest


@pytest.mark.real_data
def test_qc_report(real_data_cache, real_cfg, tmp_path):
    from cgr_gwas_qc.workflow.scripts import qc_report

    outfile = tmp_path / "qc_report.md"

    qc_report.main(
        config=real_cfg.config,
        sample_sheet_csv=real_data_cache / "dev_outputs/cgr_sample_sheet.csv",
        snp_qc_csv=real_data_cache / "dev_outputs/sample_level/snp_qc.csv",
        sample_qc_csv=real_data_cache / "dev_outputs/sample_level/sample_qc.csv",
        subject_qc_csv=real_data_cache / "dev_outputs/subject_level/subject_qc.csv",
        population_qc_csv=real_data_cache / "dev_outputs/subject_level/population_qc.csv",
        control_replicates_csv=real_data_cache
        / "dev_outputs/sample_level/concordance/InternalQcKnown.csv",
        study_replicates_csv=real_data_cache
        / "dev_outputs/sample_level/concordance/StudySampleKnown.csv",
        unexpected_replicates_csv=real_data_cache
        / "dev_outputs/sample_level/concordance/UnknownReplicates.csv",
        call_rate_png=real_data_cache / "dev_outputs/sample_level/call_rate.png",
        chrx_inbreeding_png=real_data_cache / "dev_outputs/sample_level/chrx_inbreeding.png",
        ancestry_png=real_data_cache / "dev_outputs/sample_level/ancestry.png",
        autosomal_heterozygosity_png_dir=real_data_cache
        / "dev_outputs/subject_level/autosomal_heterozygosity_plots",
        pca_png_dir=real_data_cache / "dev_outputs/subject_level/pca_plots",
        hwe_png_dir=real_data_cache / "dev_outputs/subject_level/hwe_plots",
        outfile=outfile,
    )
