import pandas as pd
import pytest
from pandas.testing import assert_frame_equal, assert_series_equal

from cgr_gwas_qc.testing import comparison, run_snakemake
from cgr_gwas_qc.testing.data import RealData


################################################################################
# Legacy Regression Tests
################################################################################
@pytest.mark.xfail(reason="I don't have a BED parser to do the comparisons.")
@pytest.mark.parametrize(
    "legacy,dev",
    [
        pytest.param(
            "legacy_outputs/subject_level/samples.bed",
            "dev_outputs/subject_level/samples.bed",
            id="legacy-samples-bed",
        ),
        pytest.param(
            "legacy_outputs/subject_level/subjects.bed",
            "dev_outputs/subject_level/subjects.bed",
            id="legacy-subjects-bed",
        ),
    ],
)
@pytest.mark.real_data
@pytest.mark.regression
def test_legacy_bed(real_data_cache, legacy, dev):
    # TODO: Implement a BED comparisons.
    comparison.assert_plink_bed_equal(real_data_cache / legacy, real_data_cache / dev)


@pytest.mark.parametrize(
    "legacy,dev",
    [
        pytest.param(
            "legacy_outputs/subject_level/samples.bim",
            "dev_outputs/subject_level/samples.bim",
            id="legacy-samples-bim",
        ),
        pytest.param(
            "legacy_outputs/subject_level/subjects.bim",
            "dev_outputs/subject_level/subjects.bim",
            id="legacy-subjects-bim",
        ),
    ],
)
@pytest.mark.real_data
@pytest.mark.regression
def test_legacy_bim(real_data_cache, legacy, dev):
    comparison.assert_plink_bim_equal(real_data_cache / legacy, real_data_cache / dev)


@pytest.mark.parametrize(
    "legacy,dev",
    [
        pytest.param(
            "legacy_outputs/subject_level/samples.fam",
            "dev_outputs/subject_level/samples.fam",
            id="legacy-samples-fam",
        ),
        pytest.param(
            "legacy_outputs/subject_level/subjects.fam",
            "dev_outputs/subject_level/subjects.fam",
            id="legacy-subjects-fam",
        ),
    ],
)
@pytest.mark.real_data
@pytest.mark.regression
def test_legacy_fam(real_data_cache, legacy, dev):
    comparison.assert_plink_fam_equal(real_data_cache / legacy, real_data_cache / dev)


@pytest.mark.regression
@pytest.mark.real_data
def test_subject_representative(real_data_cache):
    legacy = (
        pd.read_csv(real_data_cache / "legacy_outputs/subject_level/SampleUsedforSubject.csv")
        .dropna(subset=["Sample_ID"])
        .set_index("Sample_ID")
        .squeeze()
        .sort_index()
    )
    dev = (
        pd.read_csv(real_data_cache / "dev_outputs/subject_level/subject_qc.csv")
        .dropna(subset=["Sample_ID"])
        .set_index("Sample_ID")
        .Group_By_Subject_ID.sort_index()
    )

    assert_series_equal(legacy, dev, check_names=False)


@pytest.mark.regression
@pytest.mark.real_data
def test_european_subjects_mostly_overlap(real_data_cache):
    """Compare European subject lists.

    Ancestry estimation is using completely different methods, but the overall
    outcome should mostly be the same. Here I compare subject lists for European
    and ensure that there is >=95% overlap.
    """
    legacy = set(
        pd.read_csv(
            real_data_cache / "legacy_outputs/split_by_pop/EUR.keep.txt",
            delim_whitespace=True,
            header=None,
        )
        .iloc[:, 0]
        .tolist()
    )
    dev = set(
        pd.read_csv(
            real_data_cache / "dev_outputs/subject_level/European/subjects.txt",
            delim_whitespace=True,
            header=None,
        )
        .iloc[:, 0]
        .tolist()
    )

    prop_overlapping_samples = len(legacy.intersection(dev)) / len(legacy.union(dev))
    assert prop_overlapping_samples >= 0.95


@pytest.mark.regression
@pytest.mark.real_data
def test_european_autosomal_heterozygosity(real_data_cache):
    from cgr_gwas_qc.parsers import plink

    legacy = plink.read_het(
        real_data_cache / "legacy_outputs/autosomal_heterozygosity/EUR_subjects_qc.het"
    )
    dev = plink.read_het(real_data_cache / "dev_outputs/subject_level/European/subjects.het")

    overlap = legacy.index[legacy.index.isin(dev.index)]
    legacy2 = legacy.reindex(overlap)
    dev2 = dev.reindex(overlap)

    # I don't know why O(HOM) would be different
    assert_series_equal(legacy2["O_HOM"], dev2["O_HOM"], rtol=0.002)

    # E(HOM) is based on MAFs, so I do expect some differences
    assert_series_equal(legacy2["E_HOM"], dev2["E_HOM"], rtol=0.002)

    # I don't know why N(NM) would be different
    assert_series_equal(legacy2["N_NM"], dev2["N_NM"], rtol=0.002)

    # F is based on the above values so I do expect some small magnitude differences
    assert_series_equal(legacy2["F"], dev2["F"], atol=0.001)


@pytest.mark.xfail(reason="You cannot really compare PCA computed with different sets.")
@pytest.mark.regression
@pytest.mark.real_data
def test_european_pca(real_data_cache):
    from cgr_gwas_qc.parsers import eigensoft

    legacy = eigensoft.Eigenvec(
        real_data_cache / "legacy_outputs/pca/EUR_subjects.eigenvec"
    ).components
    dev = eigensoft.Eigenvec(
        real_data_cache
        / "dev_outputs/subject_level/European/subjects_unrelated_maf0.2_ld0.1.eigenvec"
    ).components

    # Pull out overlapping samples and take absolute value of all components
    # b/c sign does not mean anything in PCA.
    overlap = legacy.index[legacy.index.isin(dev.index)]
    legacy2 = legacy.reindex(overlap).applymap(abs)
    dev2 = dev.reindex(overlap).applymap(abs)

    # I don't know why O(HOM) would be different
    assert_frame_equal(legacy2, dev2, atol=0.01)


@pytest.mark.slow
@pytest.mark.regression
@pytest.mark.real_data
def test_european_hwe(real_data_cache):
    from cgr_gwas_qc.parsers import plink

    legacy = plink.read_hwe(real_data_cache / "legacy_outputs/HWP/EUR_subjects_qc.hwe")
    dev = plink.read_hwe(
        real_data_cache
        / "dev_outputs/subject_level/European/controls_unrelated_maf0.05_snps_autosome.hwe"
    )

    # Check overlap
    legacy_snp_ids = set(legacy.index)
    dev_snp_ids = set(dev.index)
    overlap = legacy_snp_ids.intersection(dev_snp_ids)
    prop_overlap = len(overlap) / len(legacy_snp_ids.union(dev_snp_ids))
    assert prop_overlap >= 0.95

    legacy2 = legacy.reindex(overlap)
    dev2 = dev.reindex(overlap)

    # These should always be the same
    assert_series_equal(legacy2["CHR"], dev2["CHR"])
    assert_series_equal(legacy2["TEST"], dev2["TEST"])

    # These values are estimated from the data, so could be very different.
    assert_series_equal(legacy2["O_HET"], dev2["O_HET"], atol=0.1)
    assert_series_equal(legacy2["E_HET"], dev2["E_HET"], atol=0.1)
    assert_series_equal(legacy2["P"], dev2["P"], atol=1)  # P is very different for a small number

    # More complicated comparisons
    # Ignore allele order
    assert all(
        ((legacy2.A1 == dev2.A1) & (legacy2.A2 == dev2.A2))
        | ((legacy2.A1 == dev2.A2) & (legacy2.A2 == dev2.A1))
    )

    # Genotypes
    legacy_geno = pd.DataFrame(
        legacy2.GENO.str.split("/").values.tolist(), columns=["A1", "het", "A2"]
    ).astype(int)
    dev_geno = pd.DataFrame(
        dev2.GENO.str.split("/").values.tolist(), columns=["A1", "het", "A2"]
    ).astype(int)
    assert_frame_equal(legacy_geno, dev_geno, atol=5)


################################################################################
# Workflow Tests
################################################################################
# -------------------------------------------------------------------------------
# Defaults
# -------------------------------------------------------------------------------
@pytest.mark.real_data
@pytest.mark.workflow
@pytest.fixture(scope="module")
def subject_qc_workflow(pytestconfig, tmp_path_factory, conda_envs):
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    tmp_path = tmp_path_factory.mktemp("subject_qc_workflow")
    conda_envs.copy_env("plink2", tmp_path)
    conda_envs.copy_env("eigensoft", tmp_path)

    (
        RealData(tmp_path)
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            module subject_qc:
                snakefile: cfg.subworkflow("subject_qc")
                config: {}

            use rule * from subject_qc
            """
        )
        .make_cgr_sample_sheet()
        .copy("dev_outputs/config.yml", "config.yml")
        .copy("dev_outputs/sample_level/sample_qc.csv", "sample_level/sample_qc.csv")
        .copy(
            "dev_outputs/sample_level/concordance/summary.csv",
            "sample_level/concordance/summary.csv",
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
    )

    run_snakemake(tmp_path, keep_temp=True)

    return tmp_path


@pytest.mark.parametrize(
    "filename",
    [
        pytest.param("subject_level/samples.bed", id="samples-bed"),
        pytest.param("subject_level/samples.bim", id="samples-bim"),
        pytest.param("subject_level/samples.fam", id="samples-fam"),
        pytest.param("subject_level/subjects.bed", id="subject-bed"),
        pytest.param("subject_level/subjects.bim", id="subject-bim"),
        pytest.param("subject_level/subjects.fam", id="subject-fam"),
        pytest.param("subject_level/concordance.csv", id="concordance"),
        pytest.param("subject_level/population_qc.csv", id="population_qc"),
    ],
)
@pytest.mark.real_data
@pytest.mark.workflow
def test_subject_qc_file_hashes(real_data_cache, subject_qc_workflow, filename):
    assert comparison.file_hashes_equal(
        real_data_cache / "dev_outputs" / filename, subject_qc_workflow / filename,
    )
