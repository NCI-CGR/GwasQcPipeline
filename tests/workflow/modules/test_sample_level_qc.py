import shutil

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from cgr_gwas_qc import load_config
from cgr_gwas_qc.testing import chdir, file_hashes_equal, make_snakefile, run_snakemake
from cgr_gwas_qc.testing.data import RealData


################################################################################
# Additional Filters and Conversions
################################################################################
@pytest.mark.real_data
@pytest.mark.workflow
@pytest.fixture(scope="session")
def graf_inputs(tmp_path_factory, conda_envs):
    # GIVEN: samples after call rate 2 and GRAF
    tmp_path = tmp_path_factory.mktemp("graf_inputs")
    conda_envs.copy_env("graf", tmp_path)
    (
        RealData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples.bed",
            "sample_level/call_rate_2/samples.bed",
        )
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples.bim",
            "sample_level/call_rate_2/samples.bim",
        )
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples.fam",
            "sample_level/call_rate_2/samples.fam",
        )
        .make_config(workflow_params={"subject_id_to_use": "PI_Subject_ID"})
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("sample_level_qc.smk")

            rule all:
                input:
                    "sample_level/call_rate_2/samples_1kg_rsID.fpg"
            """
        )
    )

    # WHEN: I run snakemake to create the GRAF FPG
    run_snakemake(tmp_path)

    return tmp_path


@pytest.mark.real_data
@pytest.mark.workflow
def test_update_snps_to_1kg_rsID(graf_inputs):
    # GIVEN: outputs from update_snps_to_1kg_rsID

    # THEN: The bim file should exist
    assert (graf_inputs / "sample_level/call_rate_2/samples_1kg_rsID.bim").exists()

    # and be different from the original BIM
    assert not file_hashes_equal(
        graf_inputs / "sample_level/call_rate_2/samples_1kg_rsID.bim",
        graf_inputs / "sample_level/call_rate_2/samples.bim",
    )

    # The BED and FAM files should be identical to the original file
    assert file_hashes_equal(
        graf_inputs / "sample_level/call_rate_2/samples_1kg_rsID.bed",
        graf_inputs / "sample_level/call_rate_2/samples.bed",
    )
    assert file_hashes_equal(
        graf_inputs / "sample_level/call_rate_2/samples_1kg_rsID.fam",
        graf_inputs / "sample_level/call_rate_2/samples.fam",
    )


@pytest.mark.real_data
@pytest.mark.workflow
def test_graf_extract_fingerprint_snps(graf_inputs):
    # GIVEN: the log from graf_extract_finger_print_snps
    graf_log = (graf_inputs / "sample_level/call_rate_2/samples_1kg_rsID.fpg.log").read_text()

    # THEN: these values should be in the log
    assert "Found total 206 samples" in graf_log
    assert (
        "3386 SNPs are FP SNPs." in graf_log
    )  # if I don't update 1kg rsIDs then there are only 2611
    assert (
        "0 SNPs need to be flipped" in graf_log
    )  # if I don't update 1kg rsIDs then there are 1327


################################################################################
# Sample Contamination
################################################################################


################################################################################
# Sample/Replicate Concordance
################################################################################
@pytest.mark.real_data
@pytest.mark.workflow
@pytest.mark.regression
def test_sample_concordance_plink(tmp_path):
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


@pytest.mark.real_data
@pytest.mark.workflow
@pytest.mark.regression
def test_sample_concordance_graf(tmp_path, graf_inputs, conda_envs):
    conda_envs.copy_env("graf", tmp_path)
    shutil.copytree(graf_inputs, tmp_path, dirs_exist_ok=True)
    make_snakefile(
        tmp_path,
        """
        from cgr_gwas_qc import load_config

        cfg = load_config()

        include: cfg.modules("sample_level_qc.smk")

        rule all:
            input:
                "sample_level/concordance/graf_relatedness.txt"
        """,
    )

    run_snakemake(tmp_path)

    # THEN: Samples marked as ID should also be marked as dups by plink ibd
    obs_graf = pd.read_csv(
        tmp_path / "sample_level/concordance/graf_relatedness.txt", sep="\t", comment="#"
    )
    obs_dups = {
        tuple(sorted(x))
        for x in obs_graf.query("`geno relation` == 'ID'")[["sample1", "sample2"]].itertuples(
            index=False
        )
    }  # Set of tuples of Sample_IDs that are identical

    exp_plink = pd.read_csv(RealData() / "production_outputs/concordance/KnownReplicates.csv")
    exp_dups = {
        tuple(sorted(x))
        for x in exp_plink.dropna()[["Sample_ID1", "Sample_ID2"]].itertuples(index=False)
    }  # Set of tuples of Sample_IDs that are concordant

    assert obs_dups == exp_dups


################################################################################
# Sample Level Ancestry
################################################################################
@pytest.mark.real_data
@pytest.mark.workflow
@pytest.mark.regression
def test_graf_ancestry(tmp_path, graf_inputs, conda_envs):
    conda_envs.copy_env("graf", tmp_path)
    shutil.copytree(graf_inputs, tmp_path, dirs_exist_ok=True)
    make_snakefile(
        tmp_path,
        """
        from cgr_gwas_qc import load_config

        cfg = load_config()

        include: cfg.modules("sample_level_qc.smk")

        rule all:
            input:
                "sample_level/ancestry/graf_ancestry.txt"
        """,
    )

    run_snakemake(tmp_path)

    # THEN: GRAF ancestry percentages should be close to the legacy SNPweights
    # ancestry percentages
    obs_ = (
        pd.read_csv(tmp_path / "sample_level/ancestry/graf_ancestry.txt", sep="\t")
        .rename(
            {
                "Subject": "Sample_ID",
                "P_f (%)": "pct_AFR",
                "P_e (%)": "pct_EUR",
                "P_a (%)": "pct_ASN",
            },
            axis=1,
        )
        .set_index("Sample_ID")
        .reindex(["pct_AFR", "pct_EUR", "pct_ASN"], axis=1)
        .sort_index()
    )

    exp_ = (
        pd.read_csv(RealData() / "production_outputs/snpweights/samples.snpweights.csv")
        .rename({"ID": "Sample_ID"}, axis=1)
        .set_index("Sample_ID")
        .assign(
            pct_AFR=lambda x: x.AFR * 100,
            pct_EUR=lambda x: x.EUR * 100,
            pct_ASN=lambda x: x.ASN * 100,
        )
        .reindex(["pct_AFR", "pct_EUR", "pct_ASN"], axis=1)
        .sort_index()
    )
    df = obs_.join(exp_, lsuffix="_obs", rsuffix="_exp")
    corr_ = df.corr()

    assert corr_.loc["pct_EUR_obs", "pct_EUR_exp"] >= 0.95
    assert corr_.loc["pct_AFR_obs", "pct_AFR_exp"] >= 0.95
    assert corr_.loc["pct_ASN_obs", "pct_ASN_exp"] >= 0.95
