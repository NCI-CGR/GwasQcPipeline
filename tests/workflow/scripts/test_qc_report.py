import shutil

import pandas as pd
import pytest

from cgr_gwas_qc.testing.data import RealData


@pytest.mark.real_data
@pytest.fixture
def control_replicates(tmp_path):
    outfile = tmp_path / "InternalQcKnown.csv"

    (
        pd.read_csv(RealData() / "production_outputs/concordance/InternalQcKnown.csv")
        .rename({"Concordance": "concordance"}, axis=1)
        .to_csv(outfile)
    )

    return outfile


@pytest.mark.real_data
@pytest.fixture
def study_replicates(tmp_path):
    outfile = tmp_path / "StudySampleKnown.csv"
    (
        pd.read_csv(RealData() / "production_outputs/concordance/StudySampleKnown.csv")
        .rename({"Concordance": "concordance"}, axis=1)
        .to_csv(outfile)
    )

    return outfile


@pytest.mark.real_data
@pytest.fixture
def unknown_replicates(tmp_path):
    outfile = tmp_path / "UnknownReplicates.csv"
    (
        pd.read_csv(RealData() / "production_outputs/concordance/UnknownReplicates.csv")
        .rename({"Concordance": "concordance"}, axis=1)
        .to_csv(outfile)
    )

    return outfile


@pytest.mark.real_data
def test_qc_report(
    real_cfg,
    snp_qc,
    sample_qc,
    population_qc,
    control_replicates,
    study_replicates,
    unknown_replicates,
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
        snp_qc_csv=snp_qc,
        sample_qc_csv=sample_qc,
        population_qc_csv=population_qc,
        control_replicates_csv=control_replicates,
        study_replicates_csv=study_replicates,
        unexpected_replicates_csv=unknown_replicates,
        call_rate_png=fake_image,
        chrx_inbreeding_png=fake_image,
        ancestry_png=fake_image,
        autosomal_heterozygosity_png_dir=png_dir,
        pca_png_dir=png_dir,
        hwe_png_dir=png_dir,
        outfile=outfile,
    )
