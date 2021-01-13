import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from cgr_gwas_qc import load_config
from cgr_gwas_qc.testing import chdir, run_snakemake
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.real_data
@pytest.mark.workflow
def test_sample_concordance(tmp_path):
    # GIVEN: Real data
    data_cache = (
        RealData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .make_config(workflow_params={"subject_id_to_use": "PI_Subject_ID"})
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("sample_level_qc.smk")

            rule all:
                input:
                    "sample_level/concordance/KnownReplicates.csv",
                    "sample_level/concordance/InternalQcKnown.csv",
                    "sample_level/concordance/StudySampleKnown.csv",
                    "sample_level/concordance/UnknownReplicates.csv",
            """
        )
    )

    with chdir(tmp_path):
        cfg = load_config()
        maf, ld = cfg.config.software_params.maf_for_ibd, cfg.config.software_params.ld_prune_r2

    (
        data_cache.copy(
            "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss",
            "sample_level/call_rate_2/samples.imiss",
        ).copy(
            "production_outputs/ibd/samples.genome",
            f"sample_level/call_rate_2/samples_maf{maf}_ld{ld}.genome",
        )
    )

    # WHEN: we run snakemake looking for the replicate output files.
    run_snakemake(tmp_path)

    # THEN:
    def _compare(obs_, exp_):
        assert_frame_equal(
            (
                pd.read_csv(obs_)
                .sort_values(["Subject_ID", "Sample_ID1", "Sample_ID2"])
                .reset_index(drop=True)
            ),
            (
                pd.read_csv(exp_)
                .rename({"Concordance": "concordance"}, axis=1)
                .dropna()  # The legacy workflow incorrectly call samples concordant if they have NAs
                .sort_values(["Subject_ID", "Sample_ID1", "Sample_ID2"])
                .reset_index(drop=True)
            ),
            check_exact=False,
        )

    _compare(
        tmp_path / "sample_level/concordance/KnownReplicates.csv",
        data_cache / "production_outputs/concordance/KnownReplicates.csv",
    )

    _compare(
        tmp_path / "sample_level/concordance/InternalQcKnown.csv",
        data_cache / "production_outputs/concordance/InternalQcKnown.csv",
    )

    _compare(
        tmp_path / "sample_level/concordance/StudySampleKnown.csv",
        data_cache / "production_outputs/concordance/StudySampleKnown.csv",
    )
