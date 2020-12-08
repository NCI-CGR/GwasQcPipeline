import pandas as pd
import pytest
from numpy import isclose

from cgr_gwas_qc.parsers.illumina.adpc import AdpcReader
from cgr_gwas_qc.testing import file_hashes_equal, file_rows_almost_equal, run_snakemake
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_median_idat_intensity(tmp_path, conda_envs):
    # GIVEN: real data repo and the illuminaio conda env
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
                input:
                    cfg.expand("sample_filters/median_idat_intensity/{Sample_ID}.{SentrixBarcode_A}.{SentrixPosition_A}.csv")
            """
        )
    )

    # WHEN: I run snakemake to build the median idat intensities
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: The observed and expected median values will match
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


@pytest.mark.real_data
@pytest.mark.workflow
def test_agg_median_idat_intensity(tmp_path, conda_envs):
    # GIVEN: A real data with IDATs and a Illuminaio conda env.
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
    # WHEN: I run snakemake to build an aggregated table of median IDAT intensities
    run_snakemake(tmp_path)

    # THEN:
    obs_table = pd.read_csv(tmp_path / "sample_filters/agg_median_idat_intensity.csv")

    # The output table has the following columns
    assert sorted(obs_table.columns) == sorted(["Sample_ID", "Chip_ID", "median_intensity"])

    # The output table has the same number of samples as the sample sheet
    obs_n_samples = obs_table.shape[0]
    exp_n_samples = data_store.ss.data.shape[0]
    assert obs_n_samples == exp_n_samples


@pytest.mark.xfail(reason="Very small floats are causing NaN. See issue #31")
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

    # THEN:
    for r in data_store.ss.data.itertuples(index=False):
        # The binary adpc.bin files should be the same accounting for slight
        # difference between floating point number
        obs_adpc = tmp_path / f"sample_filters/convert_gtc_to_illumina_adpc/{r.Sample_ID}.adpc.bin"
        exp_adpc = data_store / f"production_outputs/contam/{r.Sample_ID}.adpc.bin"
        with AdpcReader(obs_adpc) as obs, AdpcReader(exp_adpc) as expect:
            for exp_row, obs_row in zip(expect, obs):
                # each int should be the exact same
                assert exp_row.x_raw == obs_row.x_raw
                assert exp_row.y_raw == obs_row.y_raw
                assert exp_row.genotype == obs_row.genotype

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
    # WHEN: I run snakemake to convert GTC to ADPC.bin
    run_snakemake(tmp_path)

    # THEN:
    obs_abf = tmp_path / "sample_filters/GSAMD-24v1-0_20011747_A1.AF.abf.txt"
    exp_abf = data_store / "production_outputs/GSAMD-24v1-0_20011747_A1.AF.abf.txt"
    assert file_rows_almost_equal(obs_abf, exp_abf, fuzzy_col=1, sep="\t", header=True)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def contamination_test(tmp_path, conda_envs):
    # GIVEN:
    conda_envs.copy_env("verifyidintensity", tmp_path)
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
                    cfg.expand("sample_filters/contamination_test/{Sample_ID}.contam.out")
            """
        )
    )

    (
        data_store.copy(
            "production_outputs/contam", "sample_filters/convert_gtc_to_illumina_adpc"
        ).copy(
            "production_outputs/GSAMD-24v1-0_20011747_A1.AF.abf.txt",
            "sample_filters/GSAMD-24v1-0_20011747_A1.AF.abf.txt",
        )
    )

    # WHEN:
    run_snakemake(tmp_path, keep_temp=True)

    # THEN:
    for r in data_store.ss.data.itertuples(index=False):
        obs_file = tmp_path / f"sample_filters/contamination_test/{r.Sample_ID}.contam.out"
        exp_file = (
            data_store / f"production_outputs/one_samp_b_1000g_contam/{r.Sample_ID}.contam.out"
        )
        assert file_hashes_equal(obs_file, exp_file)
