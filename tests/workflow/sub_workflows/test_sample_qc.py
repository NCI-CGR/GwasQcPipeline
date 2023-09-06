import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from cgr_gwas_qc.testing import comparison, run_snakemake
from cgr_gwas_qc.testing.data import RealData


################################################################################
# Legacy Regression Tests
################################################################################
@pytest.mark.skip(reason="I don't have a BED parser to do the comparisons.")
@pytest.mark.parametrize(
    "legacy,dev",
    [
        pytest.param(
            "legacy_outputs/plink_filter_call_rate_1/samples.bed",
            "dev_outputs/sample_level/call_rate_1/samples.bed",
            id="cr1-bed",
        ),
        pytest.param(
            "legacy_outputs/plink_filter_call_rate_2/samples.bed",
            "dev_outputs/sample_level/call_rate_2/samples.bed",
            id="cr2-bed",
        ),
        pytest.param(
            "legacy_outputs/ld_prune/samples.bed",
            "dev_outputs/sample_level/call_rate_2/samples_maf0.2_ld0.1_pruned.bed",
            id="ld_pruned-bed",
        ),
    ],
)
@pytest.mark.real_data
@pytest.mark.regression
def test_legacy_samples_bed(real_data_cache, legacy, dev):
    """sample_level/samples.bed"""
    # TODO: Implement a BED comparisons.
    comparison.assert_plink_bed_equal(real_data_cache / legacy, real_data_cache / dev)


@pytest.mark.parametrize(
    "legacy,dev",
    [
        pytest.param(
            "legacy_outputs/plink_filter_call_rate_1/samples.bim",
            "dev_outputs/sample_level/call_rate_1/samples.bim",
            id="cr1-bim",
        ),
        pytest.param(
            "legacy_outputs/plink_filter_call_rate_2/samples.bim",
            "dev_outputs/sample_level/call_rate_2/samples.bim",
            id="cr2-bim",
        ),
        pytest.param(
            "legacy_outputs/ld_prune/samples.bim",
            "dev_outputs/sample_level/call_rate_2/samples_maf0.2_ld0.1_pruned.bim",
            id="ld_pruned-bim",
        ),
    ],
)
@pytest.mark.real_data
@pytest.mark.regression
def test_legacy_samples_bim(real_data_cache, legacy, dev):
    """sample_level/samples.bim"""
    comparison.assert_plink_bim_equal(real_data_cache / legacy, real_data_cache / dev)


@pytest.mark.parametrize(
    "legacy,dev",
    [
        pytest.param(
            "legacy_outputs/plink_filter_call_rate_1/samples.fam",
            "dev_outputs/sample_level/call_rate_1/samples.fam",
            id="cr1-fam",
        ),
        pytest.param(
            "legacy_outputs/plink_filter_call_rate_2/samples.fam",
            "dev_outputs/sample_level/call_rate_2/samples.fam",
            id="cr2-fam",
        ),
        pytest.param(
            "legacy_outputs/ld_prune/samples.fam",
            "dev_outputs/sample_level/call_rate_2/samples_maf0.2_ld0.1_pruned.fam",
            id="ld_pruned-fam",
        ),
    ],
)
@pytest.mark.real_data
@pytest.mark.regression
def test_legacy_samples_fam(real_data_cache, legacy, dev):
    """sample_level/samples.fam"""
    comparison.assert_plink_fam_equal(real_data_cache / legacy, real_data_cache / dev)


@pytest.mark.parametrize(
    "legacy,dev",
    [
        pytest.param(
            "legacy_outputs/concordance/KnownReplicates.csv",
            "dev_outputs/sample_level/concordance/KnownReplicates.csv",
            id="KnownReplicates",
        ),
        pytest.param(
            "legacy_outputs/concordance/InternalQcKnown.csv",
            "dev_outputs/sample_level/concordance/InternalQcKnown.csv",
            id="InternalQcKnown",
        ),
        pytest.param(
            "legacy_outputs/concordance/StudySampleKnown.csv",
            "dev_outputs/sample_level/concordance/StudySampleKnown.csv",
            id="StudySampleKnown",
        ),
    ],
)
@pytest.mark.real_data
@pytest.mark.regression
def test_legacy_split_sample_concordance(real_data_cache, legacy, dev):
    legacy_df = pd.read_csv(real_data_cache / legacy).set_index(
        ["Subject_ID", "Sample_ID1", "Sample_ID2"]
    )

    dev_df = (
        pd.read_csv(real_data_cache / dev)
        .rename({"PLINK_PI_HAT": "PI_HAT", "PLINK_concordance": "Concordance"}, axis=1)
        .set_index(["Subject_ID", "Sample_ID1", "Sample_ID2"])
        .reindex(["PI_HAT", "Concordance"], axis=1)
    )

    assert_frame_equal(legacy_df, dev_df, check_exact=False, check_like=True)


@pytest.mark.real_data
@pytest.mark.regression
def test_legacy_ancestry(real_data_cache):
    legacy = (
        pd.read_csv(real_data_cache / "legacy_outputs/snpweights/samples.snpweights.csv")
        .rename({"ID": "Sample_ID"}, axis=1)
        .set_index("Sample_ID")
        .assign(
            pct_AFR=lambda x: x.AFR * 100,
            pct_EUR=lambda x: x.EUR * 100,
            pct_ASN=lambda x: x.ASN * 100,
        )
        .reindex(["pct_AFR", "pct_EUR", "pct_ASN"], axis=1)
    )

    dev = (
        pd.read_csv(
            real_data_cache / "dev_outputs/sample_level/ancestry/graf_ancestry.txt", sep="\t"
        )
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

    df = legacy.join(dev, lsuffix="_legacy", rsuffix="_dev")
    corr_ = df.corr()

    assert corr_.loc["pct_EUR_legacy", "pct_EUR_dev"] >= 0.95
    assert corr_.loc["pct_AFR_legacy", "pct_AFR_dev"] >= 0.95
    assert corr_.loc["pct_ASN_legacy", "pct_ASN_dev"] >= 0.95


@pytest.mark.real_data
@pytest.mark.regression
def test_legacy_sample_qc_table(real_data_cache):
    comparison.assert_legacy_dev_sample_qc_equal(
        real_data_cache / "legacy_outputs/all_sample_qc.csv",
        real_data_cache / "dev_outputs/sample_level/sample_qc.csv",
    )


@pytest.mark.parametrize(
    "legacy,dev",
    [
        pytest.param(
            "legacy_outputs/ibd/samples.genome",
            "dev_outputs/sample_level/call_rate_2/samples_maf0.2_ld0.1.genome",
            id="legacy_plink_ibd",
        ),
        pytest.param(
            "legacy_outputs/remove_qc_fail/LowCR.txt",
            "dev_outputs/sample_level/qc_failures/low_call_rate.txt",
            id="legacy_low_cr",
        ),
        pytest.param(
            "legacy_outputs/remove_qc_fail/contaminated.txt",
            "dev_outputs/sample_level/qc_failures/contaminated.txt",
            id="legacy_contaminated",
        ),
        pytest.param(
            "legacy_outputs/remove_qc_fail/SexDiscordant.txt",
            "dev_outputs/sample_level/qc_failures/sex_discordant.txt",
            id="legacy_sex_discord",
        ),
        pytest.param(
            "legacy_outputs/remove_qc_fail/ExpectedRepDiscordant.txt",
            "dev_outputs/sample_level/qc_failures/replicate_discordant.txt",
            id="legacy_rep_discord",
        ),
        pytest.param(
            "legacy_outputs/subject_level/internalControls.txt",
            "dev_outputs/sample_level/internal_controls.txt",
            id="legacy_internal_control",
        ),
    ],
)
@pytest.mark.real_data
@pytest.mark.regression
def test_legacy_compare_file_hashes(real_data_cache, legacy, dev):
    comparison.file_hashes_equal(real_data_cache / legacy, real_data_cache / dev)


################################################################################
# Workflow Tests
################################################################################
@pytest.mark.real_data
@pytest.mark.workflow
@pytest.fixture(scope="module")
def sample_qc_workflow(pytestconfig, tmp_path_factory, conda_envs):
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    tmp_path = tmp_path_factory.mktemp("sample_qc_workflow")
    conda_envs.copy_env("plink2", tmp_path)
    conda_envs.copy_env("graf", tmp_path)
    conda_envs.copy_env("king", tmp_path)

    (
        RealData(tmp_path)
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            module sample_qc:
                snakefile: cfg.subworkflow("sample_qc")
                config: {}

            use rule * from sample_qc
            """
        )
        .make_cgr_sample_sheet()
        .copy("dev_outputs/config.yml", tmp_path / "config.yml")
        .copy("dev_outputs/sample_level/samples.bed", "sample_level/samples.bed")
        .copy("dev_outputs/sample_level/samples.bim", "sample_level/samples.bim")
        .copy("dev_outputs/sample_level/samples.fam", "sample_level/samples.fam")
        .copy(
            "dev_outputs/sample_level/call_rate_1/samples.bed",
            "sample_level/call_rate_1/samples.bed",
        )
        .copy(
            "dev_outputs/sample_level/call_rate_1/samples.bim",
            "sample_level/call_rate_1/samples.bim",
        )
        .copy(
            "dev_outputs/sample_level/call_rate_1/samples.fam",
            "sample_level/call_rate_1/samples.fam",
        )
        .copy(
            "dev_outputs/sample_level/call_rate_2/samples.bed",
            "sample_level/call_rate_2/samples.bed",
        )
        .copy(
            "dev_outputs/sample_level/call_rate_2/samples.bim",
            "sample_level/call_rate_2/samples.bim",
        )
        .copy(
            "dev_outputs/sample_level/call_rate_2/samples.fam",
            "sample_level/call_rate_2/samples.fam",
        )
        .copy(
            "dev_outputs/sample_level/call_rate_2/samples_1kg_rsID.bed",
            "sample_level/call_rate_2/samples_1kg_rsID.bed",
        )
        .copy(
            "dev_outputs/sample_level/call_rate_2/samples_1kg_rsID.bim",
            "sample_level/call_rate_2/samples_1kg_rsID.bim",
        )
        .copy(
            "dev_outputs/sample_level/call_rate_2/samples_1kg_rsID.fam",
            "sample_level/call_rate_2/samples_1kg_rsID.fam",
        )
        .copy(
            "dev_outputs/sample_level/call_rate_2/samples_1kg_rsID.csv",
            "sample_level/call_rate_2/samples_1kg_rsID.csv",
        )
        .copy(
            "dev_outputs/sample_level/contamination/verifyIDintensity.csv",
            "sample_level/contamination/verifyIDintensity.csv",
        )
        .copy(
            "dev_outputs/sample_level/contamination/median_idat_intensity.csv",
            "sample_level/contamination/median_idat_intensity.csv",
        )
    )

    run_snakemake(tmp_path, keep_temp=True)

    return tmp_path


@pytest.mark.parametrize(
    "filename",
    [
        pytest.param("sample_level/call_rate_2/samples_maf0.2_ld0.1.genome", id="plink_ibd"),
        pytest.param("sample_level/concordance/plink.csv", id="plink_concordance"),
        pytest.param("sample_level/concordance/king.kin0", id="king_concordance"),
        pytest.param("sample_level/concordance/graf.tsv", id="graf_concordance"),
        pytest.param("sample_level/concordance/KnownReplicates.csv", id="KnownReplicates"),
        pytest.param("sample_level/concordance/InternalQcKnown.csv", id="InternalQcKnown"),
        pytest.param("sample_level/concordance/StudySampleKnown.csv", id="StudySampleKnown"),
        pytest.param("sample_level/concordance/UnknownReplicates.csv", id="UnknownReplicates"),
        pytest.param("sample_level/snp_qc.csv", id="snp_qc_table"),
        pytest.param("sample_level/sample_qc.csv", id="sample_qc_table"),
        pytest.param("sample_level/qc_failures/low_call_rate.txt", id="low_cr"),
        pytest.param("sample_level/qc_failures/contaminated.txt", id="contaminated"),
        pytest.param("sample_level/qc_failures/sex_discordant.txt", id="sex_discord"),
        pytest.param("sample_level/qc_failures/replicate_discordant.txt", id="rep_discord"),
        pytest.param("sample_level/internal_controls.txt", id="internal_control"),
        pytest.param("sample_level/call_rate.png", id="call_rate_png"),
        pytest.param("sample_level/ancestry.png", id="ancestry_png"),
        pytest.param(
            "sample_level/chrx_inbreeding.png", id="ancestry_png", marks=pytest.mark.xfail
        ),
    ],
)
@pytest.mark.real_data
@pytest.mark.workflow
def test_sample_qc_compare_hashes(real_data_cache, sample_qc_workflow, filename):
    assert comparison.file_hashes_equal(
        real_data_cache / "dev_outputs" / filename, sample_qc_workflow / filename,
    )


@pytest.mark.real_data
@pytest.mark.workflow
def test_graf_ancestry(real_data_cache, sample_qc_workflow):
    """Table order can change"""
    dev = pd.read_csv(
        real_data_cache / "dev_outputs/sample_level/ancestry/graf_ancestry.txt",
        sep="\t",
        index_col=0,
    )
    snake = pd.read_csv(
        sample_qc_workflow / "sample_level/ancestry/graf_ancestry.txt", sep="\t", index_col=0
    )

    assert_frame_equal(dev, snake, check_like=True)


@pytest.mark.real_data
@pytest.mark.workflow
def test_summary_stats(real_data_cache, sample_qc_workflow):
    """Remove paths before comparison."""
    dev = (real_data_cache / "dev_outputs/sample_level/summary_stats.txt").read_text()
    snake = (sample_qc_workflow / "sample_level/summary_stats.txt").read_text()

    res = [
        dev_row == snake_row
        for dev_row, snake_row in zip(dev.splitlines(), snake.splitlines())
        if "/" not in dev_row  # skip rows with paths
    ]

    assert all(res)
