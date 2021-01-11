import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from cgr_gwas_qc.testing import run_snakemake
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.real_data
@pytest.mark.workflow
def test_subset_all_qc(tmp_path):
    # GIVEN: The `config.yml`, `Snakefile`, and `sample_qc_report/all_samples.csv`
    data_cache = (
        RealData(tmp_path)
        .add_sample_sheet()
        .copy("production_outputs/all_sample_qc.csv", "sample_qc_report/all_samples.csv")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("docs.smk")

            rule all:
                input:
                    "word_doc/all_sample_qc.csv"
            """
        )
    )

    # WHEN: I run snakemake to build `word_doc/all_sample_qc.csv`
    run_snakemake(tmp_path)

    # THEN: The outputs match the expected outputs from the legacy workflow.
    obs_ = pd.read_csv(tmp_path / "word_doc/all_sample_qc.csv")
    exp_ = pd.read_csv(data_cache / "production_outputs/word_doc/all_sample_qc.csv")
    assert_frame_equal(obs_, exp_)
