import shutil

import pytest


@pytest.mark.real_data
def test_qc_report(
    real_cfg,
    snp_qc_csv,
    sample_qc_csv,
    population_qc_csv,
    split_sample_concordance_tables,
    fake_image,
    tmp_path,
):
    from cgr_gwas_qc.workflow.scripts import qc_report

    outfile = tmp_path / "qc_report.md"

    # Set up some place holder PNGs for PCA, AutoHet, and HWE
    png_dir = tmp_path / "fake_png_dir"
    png_dir.mkdir()
    shutil.copyfile(fake_image, png_dir / "EUR.png")
    shutil.copyfile(fake_image, png_dir / "AFR.png")

    qc_report.main(
        config=real_cfg.config,
        sample_sheet_csv=real_cfg.root / "cgr_sample_sheet.csv",
        snp_qc_csv=snp_qc_csv,
        sample_qc_csv=sample_qc_csv,
        population_qc_csv=population_qc_csv,
        control_replicates_csv=split_sample_concordance_tables / "InternalQcKnown.csv",
        study_replicates_csv=split_sample_concordance_tables / "StudySampleKnown.csv",
        unexpected_replicates_csv=split_sample_concordance_tables / "UnknownReplicates.csv",
        call_rate_png=fake_image,
        chrx_inbreeding_png=fake_image,
        ancestry_png=fake_image,
        autosomal_heterozygosity_png_dir=png_dir,
        pca_png_dir=png_dir,
        hwe_png_dir=png_dir,
        outfile=outfile,
    )
