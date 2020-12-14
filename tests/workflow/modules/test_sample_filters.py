import re
from typing import Tuple

import numpy as np
import pandas as pd
import pytest
from numpy import isclose

from cgr_gwas_qc.parsers.illumina.adpc import AdpcReader
from cgr_gwas_qc.testing import file_hashes_equal, file_rows_almost_equal, run_snakemake
from cgr_gwas_qc.testing.data import FakeData, RealData


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_agg_median_idat_intensity(tmp_path, conda_envs):
    """Runs the two rules to process IDAT intensities.

        - median_idat_intensity
        - agg_median_idat_intensity

    I cannot test `agg_median_idat_intensity` independently. Instead of
    running `median_idat_intensity` in two different tests, I decided to test
    both rules here.
    """
    # GIVEN: Real data with per sample IDATs and a conda env with Illuminaio
    conda_envs.copy_env("illuminaio", tmp_path)
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("common.smk")
            include: cfg.rules("sample_filters.smk")

            rule all:
                input: "sample_filters/agg_median_idat_intensity.csv"
            """
        )
    )

    # WHEN: I run snakemake to calculate median IDAT intensity and aggregate these results
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: The observed and expected median intensity values values will be
    # identical for each sample
    for r in data_store.ss.data.itertuples(index=False):
        obs_file = (
            tmp_path
            / f"sample_filters/median_idat_intensity/{r.Sample_ID}.{r.SentrixBarcode_A}.{r.SentrixPosition_A}.csv"
        )
        obs_counts = pd.read_csv(obs_file).median_intensity[0]

        exp_file = (
            data_store
            / f"production_outputs/idat_intensity/{r.SentrixBarcode_A}_{r.SentrixPosition_A}.intensity.txt"
        )
        exp_counts = int(exp_file.read_text().strip())

        assert obs_counts == exp_counts

    # The aggregated table will have the columns: Sample_ID, Chip_ID, and median_intensity
    obs_table = pd.read_csv(tmp_path / "sample_filters/agg_median_idat_intensity.csv")
    assert sorted(obs_table.columns) == sorted(["Sample_ID", "Chip_ID", "median_intensity"])

    # The aggregated table will have the same number of samples as the sample sheet
    obs_n_samples = obs_table.shape[0]
    exp_n_samples = data_store.ss.data.shape[0]
    assert obs_n_samples == exp_n_samples


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_convert_gtc_to_illumina_adpc(tmp_path):
    # GIVEN: Real data using GTC entry point.
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("common.smk")
            include: cfg.rules("sample_filters.smk")

            rule all:
                input:
                    cfg.expand("sample_filters/convert_gtc_to_illumina_adpc/{Sample_ID}.adpc.bin")
                    + cfg.expand("sample_filters/convert_gtc_to_illumina_adpc/{Sample_ID}.adpc.bin.numSnps.txt")
            """
        )
    )
    # WHEN: I run snakemake to convert GTC to ADPC.bin
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: The binary adpc.bin files should be the same as production. There
    # may be slight difference between floating point number
    for r in data_store.ss.data.itertuples(index=False):
        obs_adpc = tmp_path / f"sample_filters/convert_gtc_to_illumina_adpc/{r.Sample_ID}.adpc.bin"
        exp_adpc = data_store / f"production_outputs/contam/{r.Sample_ID}.adpc.bin"
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
                    assert isclose(exp_row.x_norm, obs_row.x_norm)
                    assert isclose(exp_row.y_norm, obs_row.y_norm)
                    assert isclose(exp_row.genotype_score, obs_row.genotype_score)

        obs_counts = (
            tmp_path
            / f"sample_filters/convert_gtc_to_illumina_adpc/{r.Sample_ID}.adpc.bin.numSnps.txt"
        )
        exp_counts = data_store / f"production_outputs/contam/{r.Sample_ID}.adpc.bin.numSnps.txt"
        assert file_hashes_equal(obs_counts, exp_counts)


@pytest.mark.xfail(reason="Bug causing difference. See issue #30")
@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_pull_1KG_allele_b_freq(tmp_path):
    # GIVEN: Real data using GTC entry point.
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("common.smk")
            include: cfg.rules("sample_filters.smk")

            rule all:
                input:
                    "sample_filters/{}.{}.abf.txt".format(
                        cfg.config.reference_files.illumina_manifest_file.stem,
                        cfg.config.software_params.contam_population,
                    )
            """
        )
    )
    # WHEN: I run snakemake to pull out B allele frequencies from the 1KG
    run_snakemake(tmp_path)

    # THEN: The observed and expected abf files should be the same after
    # accounting for small difference in floating point numbers.
    obs_abf = tmp_path / "sample_filters/GSAMD-24v1-0_20011747_A1.AF.abf.txt"
    exp_abf = data_store / "production_outputs/GSAMD-24v1-0_20011747_A1.AF.abf.txt"
    assert file_rows_almost_equal(obs_abf, exp_abf, fuzzy_col=1, sep="\t", header=True)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_contamination_test(tmp_path, conda_envs):
    # GIVEN: Real data with the per sample ADPC.bin files, the 1KG B allele
    # frequencies, and verifyIDintensity conda environment.
    conda_envs.copy_env("verifyidintensity", tmp_path)
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .copy("production_outputs/contam", "sample_filters/convert_gtc_to_illumina_adpc")
        .copy(
            "production_outputs/GSAMD-24v1-0_20011747_A1.AF.abf.txt",
            "sample_filters/GSAMD-24v1-0_20011747_A1.AF.abf.txt",
        )
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("common.smk")
            include: cfg.rules("sample_filters.smk")

            rule all:
                input:
                    cfg.expand("sample_filters/contamination_test/{Sample_ID}.contam.out")
            """
        )
    )

    # WHEN: I run snakemake to create per sample contamination estimates
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: The observed file should be identical to the expected file for each
    # sample.
    for r in data_store.ss.data.itertuples(index=False):
        obs_file = tmp_path / f"sample_filters/contamination_test/{r.Sample_ID}.contam.out"
        exp_file = (
            data_store / f"production_outputs/one_samp_b_1000g_contam/{r.Sample_ID}.contam.out"
        )
        assert file_hashes_equal(obs_file, exp_file)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_contamination_test_with_missing_data(tmp_path, conda_envs):
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
            "sample_filters/GSAMD-24v1-0_20011747_A1.AF.abf.txt",
        )
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("common.smk")
            include: cfg.rules("sample_filters.smk")

            rule all:
                input:
                    cfg.expand("sample_filters/contamination_test/{Sample_ID}.contam.out")
            """
        )
    )

    # WHEN: I run snakemake to create per sample contamination estimates
    # (adpc.bin) and contamination estimates (contam.out)
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: I created the adbc.bin files
    for r in data_store.ss.data.itertuples(index=False):
        assert (
            tmp_path / f"sample_filters/convert_gtc_to_illumina_adpc/{r.Sample_ID}.adpc.bin"
        ).exists()

    # The observed contamination file is identical to the expected file for each sample.
    for r in data_store.ss.data.itertuples(index=False):
        obs_file = tmp_path / f"sample_filters/contamination_test/{r.Sample_ID}.contam.out"
        exp_file = (
            data_store / f"production_outputs/one_samp_b_1000g_contam/{r.Sample_ID}.contam.out"
        )
        assert file_hashes_equal(obs_file, exp_file)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_agg_contamination_test(tmp_path, conda_envs):
    # GIVEN: Real data with per sample contamination estimates, final call
    # rates, and the Illuminaio conda env. I need Illuminaio because I have to
    # run `agg_median_idat_intensity` which does not have an equivalent in the
    # production workflow.
    conda_envs.copy_env("illuminaio", tmp_path)
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .copy("production_outputs/one_samp_b_1000g_contam", "sample_filters/contamination_test")
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss",
            "plink_filter_call_rate_2/samples.imiss",
        )
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("common.smk")
            include: cfg.rules("sample_filters.smk")

            rule all:
                input:
                    "sample_filters/agg_contamination_test.csv"
            """
        )
    )

    # WHEN: I run snakemake to build the aggregated contamination table.
    run_snakemake(tmp_path)

    # THEN: The observed and expected tables have almost identical %Mix after
    # accounting for small floating point differences
    obs_file = tmp_path / "sample_filters/agg_contamination_test.csv"
    obs = pd.read_csv(obs_file)

    exp_file = data_store / "production_outputs/all_contam/contam.csv"
    exp_ = pd.read_csv(exp_file)

    sample_ids = obs.Sample_ID.to_list()  # noqa: F841
    assert all(isclose(exp_.query("ID == @sample_ids")["%Mix"].values, obs["%Mix"].values))


@pytest.mark.workflow
def test_Sample_IDs_above_contam_threshold(tmp_path):
    # GIVEN: Fake data with the GTC entry point
    (
        FakeData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("common.smk")
            include: cfg.rules("sample_filters.smk")

            rule all:
                input:
                    "sample_filters/contaminated_samples.txt"
            """
        )
    )

    # Add a fake agg_contamination_test.csv table with one sample with a %Mix
    # over the contamination threshold (default 0.2).
    (tmp_path / "sample_filters").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            ["T0001", 0.05, 0.0, 0.0],
            ["T0002", 0.1, 0.0, 0.0],
            ["T0003", 0.2, 0.0, 0.0],
            ["T0004", 0.25, 0.0, 0.0],
        ],
        columns=["Sample_ID", "%Mix", "LLK", "LLK0"],
    ).to_csv(tmp_path / "sample_filters/agg_contamination_test.csv", index=False)

    # WHEN: I run snakemake to create a list of contaminated samples
    run_snakemake(tmp_path)

    # THEN: I expect a single sample to be flagged as contaminated.
    obs_file = tmp_path / "sample_filters/contaminated_samples.txt"
    obs = obs_file.read_text().strip().split("\n")
    assert len(obs) == 1


def parse_plink_contamination_filtering(file_name) -> Tuple[int, int]:
    m = re.findall(
        r"\n(\d+) people.*\n--remove: (\d+) people remaining", file_name.read_text(), re.DOTALL
    )
    n_start, n_remaining = m[0]
    return int(n_start), int(n_remaining)


@pytest.mark.real_data
@pytest.mark.workflow
def test_remove_contaminated_samples(tmp_path, conda_envs):
    # GIVEN: Real data with the GTC entry point, plink files after the second
    # call rate filter, and the plink conda environment
    conda_envs.copy_env("plink2", tmp_path)
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .copy("production_outputs/plink_filter_call_rate_2", "plink_filter_call_rate_2")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("common.smk")
            include: cfg.rules("sample_filters.smk")

            rule all:
                input:
                    expand("sample_filters/not_contaminated/samples.{ext}", ext=["bed", "bim",  "fam"])
            """
        )
    )

    # and a contaminated sample list with 1 sample id.
    contam_sample = data_store.ss.data.Sample_ID[0]
    (tmp_path / "sample_filters").mkdir(parents=True, exist_ok=True)
    (tmp_path / "sample_filters/contaminated_samples.txt").write_text(
        "{} {}\n".format(contam_sample, contam_sample)
    )

    # WHEN: I run snakemake to remove contaminated samples
    run_snakemake(tmp_path)

    # THEN: I expect the files will be created
    outputs = [
        tmp_path / "sample_filters/not_contaminated/samples.bed",
        tmp_path / "sample_filters/not_contaminated/samples.bim",
        tmp_path / "sample_filters/not_contaminated/samples.fam",
    ]
    assert all(file_name.exists() for file_name in outputs)

    # and the plink log shows 1 sample was filtered
    n_start, n_remaining = parse_plink_contamination_filtering(
        tmp_path / "sample_filters/not_contaminated/samples.log"
    )
    assert (n_start - n_remaining) == 1


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_plink_ibd(tmp_path, conda_envs):
    # GIVEN:
    conda_envs.copy_env("plink2", tmp_path)
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .make_config()
        .copy("production_outputs/ld_prune", "snp_filters/ld_prune")
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("common.smk")
            include: cfg.rules("sample_filters.smk")

            rule all:
                input:
                    expand("sample_filters/ibd/samples.{ext}", ext=["genome", "nosex"])
            """
        )
    )

    # WHEN:
    run_snakemake(tmp_path)

    # THEN:
    obs_file = tmp_path / "sample_filters/ibd/samples.genome"
    exp_file = data_store / "production_outputs/ibd/samples.genome"
    assert file_hashes_equal(obs_file, exp_file)


@pytest.mark.real_data
@pytest.mark.workflow
def test_sample_concordance(tmp_path):
    # GIVEN: Real data
    (
        RealData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .make_config()
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples_filter2.imiss",
            "plink_filter_call_rate_2/samples.imiss",
        )
        .copy("production_outputs/ibd/samples.genome", "sample_filters/ibd/samples.genome")
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("common.smk")
            include: cfg.rules("sample_filters.smk")

            rule all:
                input:
                    "sample_filters/concordance/KnownReplicates.csv",
                    "sample_filters/concordance/InternalQcKnown.csv",
                    "sample_filters/concordance/StudySampleKnown.csv",
                    "sample_filters/concordance/UnknownReplicates.csv",
            """
        )
    )

    # WHEN: we run snakemake looking for the replicate output files.
    run_snakemake(tmp_path)

    # THEN: these files are created.
    # Note: regression tests are executed as part of the
    # `known_concorant_samples.py` testing
    assert (tmp_path / "sample_filters/concordance/KnownReplicates.csv").exists()
    assert (tmp_path / "sample_filters/concordance/InternalQcKnown.csv").exists()
    assert (tmp_path / "sample_filters/concordance/StudySampleKnown.csv").exists()
    assert (tmp_path / "sample_filters/concordance/UnknownReplicates.csv").exists()
