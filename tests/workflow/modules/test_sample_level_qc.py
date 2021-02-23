import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from cgr_gwas_qc import load_config
from cgr_gwas_qc.parsers.illumina import AdpcReader
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
@pytest.mark.slow
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
@pytest.mark.slow
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
@pytest.mark.real_data
@pytest.mark.workflow
@pytest.fixture(scope="session")
def median_idat_intensity(tmp_path_factory, conda_envs) -> Path:
    """Aggregated median IDAT intensity.

    I changed the format of this file, so there are no legacy versions. This
    fixture will allow re-use.
    """
    # GIVEN: Real data with per sample IDATs and a conda env with Illuminaio
    tmp_path = tmp_path_factory.mktemp("idat_intensity_outputs")
    conda_envs.copy_env("illuminaio", tmp_path)
    (
        RealData(tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("common.smk")  # Need numSNPs from here
            include: cfg.modules("sample_level_qc.smk")

            rule all:
                input: "sample_level/median_idat_intensity.csv"
            """
        )
    )

    # WHEN: I run snakemake to calculate median IDAT intensity and aggregate
    # these results.
    run_snakemake(tmp_path, keep_temp=True)

    return tmp_path / "sample_level/median_idat_intensity.csv"


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_agg_median_idat_intensity(median_idat_intensity):
    # THEN: The values in the observed summary table match the values found in
    # the legacy files.
    obs_table = (
        pd.read_csv(median_idat_intensity)  # Sample_ID Chip_ID median_intensity
        .set_index("Chip_ID")
        .median_intensity
    )
    for exp_file in (RealData() / "production_outputs/idat_intensity").glob("*.txt"):
        exp_value = int(exp_file.read_text())
        exp_chip_id = exp_file.stem.split(".")[0]
        assert obs_table[exp_chip_id] == exp_value


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.slow
def test_per_sample_gtc_to_adpc(tmp_path):
    # GIVEN: Real data using GTC entry point.
    data_cache = (
        RealData(tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("common.smk")
            include: cfg.modules("sample_level_qc.smk")

            rule all:
                input:
                    cfg.expand("sample_level/per_sample_adpc/{Sample_ID}.adpc.bin")
            """
        )
    )
    # WHEN: I run snakemake to convert GTC to ADPC.bin
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: The binary adpc.bin files should be the same as production. There
    # may be slight difference between floating point number
    for r in data_cache.ss.data.itertuples(index=False):
        obs_adpc = tmp_path / f"sample_level/per_sample_adpc//{r.Sample_ID}.adpc.bin"
        exp_adpc = data_cache / f"production_outputs/contam/{r.Sample_ID}.adpc.bin"
        with AdpcReader(obs_adpc) as obs, AdpcReader(exp_adpc) as expect:
            for exp_row, obs_row in zip(expect, obs):
                # each int should be the exact same
                assert exp_row.x_raw == obs_row.x_raw
                assert exp_row.y_raw == obs_row.y_raw
                assert exp_row.genotype == obs_row.genotype

                if np.isnan(obs_row.x_norm) | np.isnan(obs_row.y_norm):
                    # Note: the legacy workflow outputs (i.e., exp_row) will
                    # have 0.0 instead of NA for x_norm and y_norm. This
                    # difference is okay, see Issue #31 for details.

                    # If normalized values are NA then the following should be true
                    assert (obs_row.x_raw, obs_row.y_raw) == (0, 0)
                    assert obs_row.genotype_score == 0.0
                    assert obs_row.genotype == 3
                else:
                    # floats should be really close
                    assert pytest.approx(exp_row.x_norm, abs=1e-6) == obs_row.x_norm
                    assert pytest.approx(exp_row.y_norm, abs=1e-6) == obs_row.y_norm
                    assert pytest.approx(exp_row.genotype_score, abs=1e-6) == obs_row.genotype_score

        obs_counts = tmp_path / f"sample_level/per_sample_num_snps/{r.Sample_ID}.txt"
        exp_counts = data_cache / f"production_outputs/contam/{r.Sample_ID}.adpc.bin.numSnps.txt"
        assert file_hashes_equal(obs_counts, exp_counts)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_per_sample_verifyIDintensity_contamination(tmp_path, conda_envs):
    # GIVEN: Real data with the per sample ADPC.bin files, the 1KG B allele
    # frequencies, and verifyIDintensity conda environment.
    conda_envs.copy_env("verifyidintensity", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .copy("production_outputs/contam", "sample_level/per_sample_adpc")
        .copy(
            "production_outputs/GSAMD-24v1-0_20011747_A1.AF.abf.txt",
            "sample_level/GSAMD-24v1-0_20011747_A1.AF.abf.txt",
        )
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("common.smk")
            include: cfg.modules("sample_level_qc.smk")

            rule all:
                input:
                    cfg.expand("sample_level/per_sample_contamination_test/{Sample_ID}.contam.out")
            """
        )
    )

    # WHEN: I run snakemake to create per sample contamination estimates
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: The observed file should be identical to the expected file for each
    # sample.
    for r in data_cache.ss.data.itertuples(index=False):
        obs_file = tmp_path / f"sample_level/per_sample_contamination_test/{r.Sample_ID}.contam.out"
        exp_file = (
            data_cache / f"production_outputs/one_samp_b_1000g_contam/{r.Sample_ID}.contam.out"
        )
        assert file_hashes_equal(obs_file, exp_file)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.slow
def test_contamination_test_with_missing_abf_values(tmp_path, conda_envs):
    """Does verifyIDintensity run with missing abf values.

    In issue #30 we found that the legacy pipeline was setting missing values
    to 0.0 when converting BPM to abf. After looking at verifyIDintensity
    code we found that these values should be `"NA"` to be excluded from the
    contamination calculation. Here I am testing that everything runs and
    generates similar outputs.
    """
    # GIVEN: Real data using GTC entry point and the per sample adpc files.
    conda_envs.copy_env("verifyidintensity", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .copy("production_outputs/contam", "sample_level/per_sample_adpc")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("common.smk")
            include: cfg.modules("sample_level_qc.smk")

            rule all:
                input:
                    cfg.expand("sample_level/per_sample_contamination_test/{Sample_ID}.contam.out")
            """
        )
    )

    # WHEN: I run snakemake to calculate contamination
    run_snakemake(tmp_path, keep_temp=True, print_stdout=True)

    # THEN: The observed contamination file is similar to the expected file for each sample.
    for r in data_cache.ss.data.itertuples(index=False):
        obs_file = tmp_path / f"sample_level/per_sample_contamination_test/{r.Sample_ID}.contam.out"
        exp_file = (
            data_cache / f"production_outputs/one_samp_b_1000g_contam/{r.Sample_ID}.contam.out"
        )

        obs = obs_file.read_text().strip().split("\n")[-1].split()
        exp_ = exp_file.read_text().strip().split("\n")[-1].split()
        assert pytest.approx(float(exp_[1]), abs=0.01) == float(obs[1])  # within 1% of each other
        assert pytest.approx(float(exp_[2]), rel=0.1) == float(obs[2])  # within 1000
        assert pytest.approx(float(exp_[3]), rel=0.1) == float(obs[3])  # within 1000


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.slow
def test_contamination_test_with_missing_adpc_values(tmp_path, conda_envs):
    """Does verifyIDintensity work when there are missing data.

    In Issue #31 I determined that there is a small difference when building
    the adpc files. The current scripts outputs normalized values as `nan`
    when raw values are 0. The legacy outputs them as 0. I want to make sure
    that verifyIDintensity does not change regardless of `nan` or 0.
    """
    # GIVEN: Real data with the 1KG B allele frequencies and verifyIDintensity
    # conda environment.
    conda_envs.copy_env("verifyidintensity", tmp_path)
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .copy(
            "production_outputs/GSAMD-24v1-0_20011747_A1.AF.abf.txt",
            "sample_level/GSAMD-24v1-0_20011747_A1.AF.abf.txt",
        )
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("common.smk")
            include: cfg.modules("sample_level_qc.smk")

            rule all:
                input:
                    cfg.expand("sample_level/per_sample_contamination_test/{Sample_ID}.contam.out")
            """
        )
    )

    # WHEN: I run snakemake to create per sample contamination estimates
    # (adpc.bin) and contamination estimates (contam.out)
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: I created the adbc.bin files
    for r in data_store.ss.data.itertuples(index=False):
        assert (tmp_path / f"sample_level/per_sample_adpc/{r.Sample_ID}.adpc.bin").exists()

    # The observed contamination file is identical to the expected file for each sample.
    for r in data_store.ss.data.itertuples(index=False):
        obs_file = tmp_path / f"sample_level/per_sample_contamination_test/{r.Sample_ID}.contam.out"
        exp_file = (
            data_store / f"production_outputs/one_samp_b_1000g_contam/{r.Sample_ID}.contam.out"
        )
        assert file_hashes_equal(obs_file, exp_file)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_agg_contamination_test(tmp_path, median_idat_intensity):
    # GIVEN: Real data with per sample contamination estimates, final call
    # rates, and the Illuminaio conda env.
    data_cache = (
        RealData(tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .copy(
            "production_outputs/one_samp_b_1000g_contam",
            "sample_level/per_sample_contamination_test/",
        )
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss",
            "sample_level/call_rate_2/samples.imiss",
        )
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("common.smk")
            include: cfg.modules("sample_level_qc.smk")

            rule all:
                input:
                    "sample_level/contamination/verifyIDintensity_contamination.csv"
            """
        )
    )
    # This file has no equivalent in the legacy workflow so I need to copy from a fixture
    shutil.copyfile(median_idat_intensity, tmp_path / "sample_level/median_idat_intensity.csv")

    # WHEN: I run snakemake to build the aggregated contamination table.
    run_snakemake(tmp_path)

    # THEN: The observed and subset of the expected have almost identical
    # values.
    obs_ = pd.read_csv(
        tmp_path / "sample_level/contamination/verifyIDintensity_contamination.csv"
    ).set_index("Sample_ID")

    exp_ = (
        pd.read_csv(data_cache / "production_outputs/all_contam/contam.csv")
        .rename({"ID": "Sample_ID"}, axis=1)
        .set_index("Sample_ID")
        .reindex(obs_.index)
    )

    assert_frame_equal(obs_, exp_)


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

            include: cfg.modules("common.smk")
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
            f"sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.genome",
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
@pytest.mark.slow
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
@pytest.mark.slow
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
