"""Tests for the entry_points module.

Each entry point is expected to create the following files:

- sample_level/samples.bed
- sample_level/samples.bim
- sample_level/samples.ped

Entry point decisions are made by the presence of different user provided
files in the config. If these options are None (or not in the file) then a
different entry point is used.

For testing entry points we need to make sure that the entry point files are
present in the working directory and that they are referenced in the config.
"""
import pytest

from cgr_gwas_qc import load_config
from cgr_gwas_qc.testing import chdir, comparison, run_snakemake
from cgr_gwas_qc.testing.comparison import file_hashes_equal
from cgr_gwas_qc.testing.data import FakeData, RealData


################################################################################
# Legacy Regression Tests
################################################################################
@pytest.mark.skip(reason="I don't have a BED parser to do the comparisons.")
@pytest.mark.real_data
@pytest.mark.regression
def test_samples_bed(real_data_cache):
    """sample_level/samples.bed"""
    # TODO: Implement a BED comparisons.
    comparison.assert_plink_bed_equal(
        real_data_cache / "legacy_outputs/plink_start/samples.bed",
        real_data_cache / "dev_outputs/sample_level/samples.bed",
    )


@pytest.mark.real_data
@pytest.mark.regression
def test_samples_bim(real_data_cache):
    """sample_level/samples.bim"""
    comparison.assert_plink_bim_equal(
        real_data_cache / "legacy_outputs/plink_start/samples.bim",
        real_data_cache / "dev_outputs/sample_level/samples.bim",
    )


@pytest.mark.real_data
@pytest.mark.regression
def test_samples_fam(real_data_cache):
    """sample_level/samples.fam"""
    comparison.assert_plink_fam_equal(
        real_data_cache / "legacy_outputs/plink_start/samples.fam",
        real_data_cache / "dev_outputs/sample_level/samples.fam",
    )


################################################################################
# Workflow Tests
################################################################################
# -------------------------------------------------------------------------------
# GTC Entry Point
# -------------------------------------------------------------------------------
@pytest.mark.real_data
@pytest.fixture(scope="module")
def gtc_entry(pytestconfig, tmp_path_factory, conda_envs):
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    tmp_path = tmp_path_factory.mktemp("gtc_entry")
    conda_envs.copy_env("plink2", tmp_path)
    (
        RealData(tmp_path, full_sample_sheet=False)
        .make_cgr_sample_sheet()
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            module entry_points:
                snakefile: cfg.subworkflow("entry_points")
                config: {"notemp": True}

            use rule * from entry_points
            """
        )
    )

    run_snakemake(tmp_path, keep_temp=True)

    return tmp_path


@pytest.mark.real_data
@pytest.mark.workflow
@pytest.mark.regression
def test_gtc_to_ped_conversion(real_data_cache, gtc_entry):
    # GIVEN: A DataRepo with GTC files and the following snakefile
    # THEN: peds and maps should be identical to production output
    legacy_peds = sorted((real_data_cache / "legacy_outputs/per_sample_plink_files").glob("*.ped"))
    snake_peds = sorted((gtc_entry / "sample_level/per_sample_plink_files").glob("*.ped"))
    assert all(
        comparison.file_hashes_equal(legacy, snake)
        for legacy, snake in zip(legacy_peds, snake_peds)
    )

    legacy_maps = sorted((real_data_cache / "legacy_outputs/per_sample_plink_files").glob("*.map"))
    snake_maps = sorted((gtc_entry / "sample_level/per_sample_plink_files").glob("*.map"))
    assert all(
        comparison.file_hashes_equal(legacy, snake)
        for legacy, snake in zip(legacy_maps, snake_maps)
    )


@pytest.mark.real_data
@pytest.mark.workflow
def test_create_gtc_merge_list(gtc_entry):
    # THEN: The merge should have one row for each sample
    with chdir(gtc_entry):
        cfg = load_config(pytest=True)

    n_samples = cfg.config.num_samples
    merge_list = (gtc_entry / "sample_level/plink_merge_list.txt").read_text().strip().split("\n")
    assert n_samples == len(merge_list)


@pytest.mark.real_data
@pytest.mark.workflow
def test_merge_gtc_sample_peds(gtc_entry):
    # The merged samples files should exist.
    assert (gtc_entry / "sample_level/samples.bed").exists()
    assert (gtc_entry / "sample_level/samples.bim").exists()
    assert (gtc_entry / "sample_level/samples.fam").exists()
    assert (gtc_entry / "sample_level/samples.nosex").exists()


@pytest.mark.real_data
@pytest.fixture(scope="module")
def gtc_grouped_entry(pytestconfig, tmp_path_factory, conda_envs):
    if not pytestconfig.getoption("--real-data"):
        pytest.skip("No real data")

    tmp_path = tmp_path_factory.mktemp("gtc_grouped_entry")
    conda_envs.copy_env("plink2", tmp_path)
    (
        RealData(tmp_path, full_sample_sheet=False)
        .make_cgr_sample_sheet()
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            module entry_points:
                snakefile: cfg.subworkflow("entry_points")
                config: {"cluster_mode": True, "notemp": True}

            use rule * from entry_points
            """
        )
    )

    run_snakemake(tmp_path, keep_temp=True)

    return tmp_path


@pytest.mark.real_data
@pytest.mark.workflow
def test_gtc_grouped_entry(gtc_grouped_entry):
    # The merged samples files should exist.
    assert (gtc_grouped_entry / "sample_level/cgroup1/samples.bed").exists()
    assert (gtc_grouped_entry / "sample_level/cgroup1/samples.bim").exists()
    assert (gtc_grouped_entry / "sample_level/cgroup1/samples.fam").exists()
    assert (gtc_grouped_entry / "sample_level/cgroup1/samples.nosex").exists()


# -------------------------------------------------------------------------------
# PED/MAP Entry Point
# -------------------------------------------------------------------------------
@pytest.mark.workflow
@pytest.fixture(scope="module")
def ped_entry(tmp_path_factory, conda_envs):
    tmp_path = tmp_path_factory.mktemp("ped_entry")
    conda_envs.copy_env("plink2", tmp_path)
    (
        FakeData(tmp_path)
        .make_cgr_sample_sheet()
        .add_user_files(entry_point="ped")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            module entry_points:
                snakefile: cfg.subworkflow("entry_points")
                config: {}

            use rule * from entry_points
            """
        )
    )

    run_snakemake(tmp_path, keep_temp=True)

    return tmp_path


@pytest.mark.regression
@pytest.mark.workflow
def test_ped_entry(fake_data_cache, ped_entry):
    """Aggregated PED entry point.

    This entry point expects an aggregated (multi-sample) PED/MAP file. I
    only have this with fake data.
    """
    # The merged samples files should be identical to the test data I have generated.
    assert file_hashes_equal(
        fake_data_cache / "plink/samples.bed", ped_entry / "sample_level/samples.bed"
    )

    assert file_hashes_equal(
        fake_data_cache / "plink/samples.bim", ped_entry / "sample_level/samples.bim"
    )

    # Fam files had strange spacing differences so just making sure the data is the same
    def parse_fam(pth):
        return sorted([row.split() for row in pth.read_text().strip().split("\n")])

    fake_fam = parse_fam(fake_data_cache / "plink/samples.fam")
    snake_fam = parse_fam(ped_entry / "sample_level/samples.fam")
    assert fake_fam == snake_fam


# -------------------------------------------------------------------------------
# BED/BIM/FAM Entry Point
# -------------------------------------------------------------------------------
@pytest.mark.workflow
@pytest.fixture(scope="module")
def bed_entry(tmp_path_factory, conda_envs):
    tmp_path = tmp_path_factory.mktemp("bed_entry")
    conda_envs.copy_env("plink2", tmp_path)
    (
        FakeData(tmp_path)
        .make_cgr_sample_sheet()
        .add_user_files(entry_point="bed")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            module entry_points:
                snakefile: cfg.subworkflow("entry_points")
                config: {}

            use rule * from entry_points
            """
        )
    )

    run_snakemake(tmp_path, keep_temp=True)

    return tmp_path


@pytest.mark.workflow
def test_bed_entry(bed_entry):
    """Aggregated BED entry point.

    This entry point expects an aggregated (multi-sample) BED/BIM/FAM file.
    Since this is the expected output of the entry point module, these files
    are just symbolically linked to the expected location
    ``sample_level/samples.{bed,bim,fam}``.
    """
    # THEN: the files should just be symlinks
    # The merged samples files should be identical to the test data I have generated.
    assert (bed_entry / "sample_level/samples.bed").is_symlink()
    assert (bed_entry / "sample_level/samples.bim").is_symlink()
    assert (bed_entry / "sample_level/samples.fam").is_symlink()
