"""Tests for the entry_points module.

Each entry point is expected to create the following files:

- plink_start/samples.bed
- plink_start/samples.bim
- plink_start/samples.ped

Entry point decisions are made by the presence of different user provided
files in the config. If these options are None (or not in the file) then a
different entry point is used.

For testing entry points we need to make sure that the entry point files are
present in the working directory and that they are referenced in the config.
"""
import pytest

from cgr_gwas_qc.testing import file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import FakeData, RealData


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_gtc_to_ped_conversion(tmp_path):
    # GIVEN: A DataRepo with GTC files and the following snakefile
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

            include: cfg.rules("entry_points.smk")

            rule all:
                input:
                    cfg.expand("per_sample_plink_files/{Sample_ID}.ped")
            """
        )
    )
    # WHEN: I run snakemake while keeping temporary outputs
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: peds and maps should be identical to production output
    obs_peds = sorted((tmp_path / "per_sample_plink_files").glob("*.ped"))
    exp_peds = sorted((tmp_path / "production_outputs/per_sample_plink_files").glob("*.ped"))
    assert all([file_hashes_equal(o, e) for o, e in zip(obs_peds, exp_peds)])

    obs_maps = sorted((tmp_path / "per_sample_plink_files").glob("*.map"))
    exp_maps = sorted((tmp_path / "production_outputs/per_sample_plink_files").glob("*.map"))
    assert all([file_hashes_equal(o, e) for o, e in zip(obs_maps, exp_maps)])


@pytest.mark.real_data
@pytest.mark.workflow
def test_create_gtc_merge_list(tmp_path):
    # GIVEN: A RealData repo with the following Snakefile
    data_repo = (
        RealData(tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("entry_points.smk")

            rule all:
                input: "plink_start/mergeList.txt"
            """
        )
    )
    # use pre-built intermediate per_sample_plink_files
    data_repo.copy("production_outputs/per_sample_plink_files", "per_sample_plink_files")

    # WHEN: I run snakemake
    run_snakemake(tmp_path)

    # THEN: The merge should have one row for each sample
    n_samples = data_repo.ss.data.shape[0]
    merge_list = (tmp_path / "plink_start/mergeList.txt").read_text().strip().split("\n")
    assert n_samples == len(merge_list)


@pytest.mark.real_data
@pytest.mark.workflow
def test_merge_gtc_sample_peds(tmp_path, conda_envs):
    # GIVEN: The plink conda env and the following RealData repo
    conda_envs.copy_env("plink2", tmp_path)
    data_repo = (
        RealData(tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("entry_points.smk")

            rule all:
                input:
                    expand("plink_start/samples.{ext}", ext=["bed", "bim", "fam", "nosex"])
            """
        )
    )
    # use pre-built intermediate per_sample_plink_files
    data_repo.copy("production_outputs/per_sample_plink_files", "per_sample_plink_files")

    # WHEN: I run snakemake
    run_snakemake(tmp_path)

    # THEN:
    # The merged samples files should exist.
    assert (tmp_path / "plink_start/samples.bed").exists()
    assert (tmp_path / "plink_start/samples.bim").exists()
    assert (tmp_path / "plink_start/samples.fam").exists()
    assert (tmp_path / "plink_start/samples.nosex").exists()

    # The same number of samples should be reported in the log as passing
    n_samples = data_repo.ss.data.shape[0]
    log = (tmp_path / "plink_start/samples.log").read_text()
    assert f"{n_samples} people pass filters and QC" in log


@pytest.mark.regression
@pytest.mark.workflow
def test_ped_entry(tmp_path, conda_envs):
    """Aggregated PED entry point.

    This entry point expects an aggregated (multi-sample) PED/MAP file. I
    only have this with fake data.
    """
    # GIVEN: The plink conda env and the following FakeData repo
    conda_envs.copy_env("plink2", tmp_path)
    data_repo = (
        FakeData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .add_user_files(entry_point="ped")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("entry_points.smk")

            rule all:
                input:
                    expand("plink_start/samples.{ext}", ext=["bed", "bim", "fam", "nosex"])
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path)

    # THEN:
    # The merged samples files should be identical to the test data I have generated.
    obs_bed = tmp_path / "plink_start/samples.bed"
    exp_bed = data_repo / "plink/samples.bed"
    assert file_hashes_equal(obs_bed, exp_bed)

    obs_bim = tmp_path / "plink_start/samples.bim"
    exp_bim = data_repo / "plink/samples.bim"
    assert file_hashes_equal(obs_bim, exp_bim)

    # Fam files had strange spacing differences so just making sure the data is the same
    def parse_fam(pth):
        return sorted([row.split() for row in pth.read_text().strip().split("\n")])

    obs_fam = parse_fam(tmp_path / "plink_start/samples.fam")
    exp_fam = parse_fam(data_repo / "plink/samples.fam")
    assert obs_fam == exp_fam


@pytest.mark.workflow
def test_bed_entry(tmp_path):
    """Aggregated BED entry point.

    This entry point expects an aggregated (multi-sample) BED/BIM/FAM file.
    Since this is the expected output of the entry point module, these files
    are just symbolically linked to the expected location
    ``plink_start/samples.{bed,bim,fam}``.
    """
    # GIVEN: The the following FakeData repo
    (
        FakeData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .add_user_files(entry_point="bed")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("entry_points.smk")

            rule all:
                input:
                    expand("plink_start/samples.{ext}", ext=["bed", "bim", "fam"])
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path)

    # THEN: the files should just be symlinks
    # The merged samples files should be identical to the test data I have generated.
    assert (tmp_path / "plink_start/samples.bed").is_symlink()
    assert (tmp_path / "plink_start/samples.bim").is_symlink()
    assert (tmp_path / "plink_start/samples.fam").is_symlink()
