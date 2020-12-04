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
from cgr_gwas_qc.testing import file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import RealData

snakefile = """\
from cgr_gwas_qc import load_config

cfg = load_config()

include: cfg.rules("entry_points.smk")

rule all:
    input:
        expand("plink_start/samples.{ext}", ext=["bed", "bim", "fam"])
"""


def test_gtc_to_ped_conversion(gtc_data_repo):
    # GIVEN: A DataRepo with GTC files and the following snakefile
    gtc_data_repo.make_snakefile(
        """
        from cgr_gwas_qc import load_config

        cfg = load_config()

        include: cfg.rules("entry_points.smk")

        rule all:
            input:
                cfg.expand("per_sample_plink_files/{Sample_ID}.ped")
        """
    )
    # WHEN: I run snakemake while keeping temporary outputs
    working_dir = gtc_data_repo.working_dir
    run_snakemake(working_dir, keep_temp=True)

    # THEN:
    n_samples = gtc_data_repo.ss.data.shape[0]
    obs_peds = sorted((working_dir / "per_sample_plink_files").glob("*.ped"))
    obs_maps = sorted((working_dir / "per_sample_plink_files").glob("*.map"))

    # I should have a PED and MAP file for each sample
    assert len(obs_peds) == n_samples
    assert len(obs_maps) == n_samples

    # If this is real data
    if isinstance(gtc_data_repo, RealData):
        # the outputs should match the production workflow
        exp_peds = sorted(
            (gtc_data_repo / "production_outputs/per_sample_plink_files").glob("*.ped")
        )
        assert all([file_hashes_equal(o, e) for o, e in zip(obs_peds, exp_peds)])

        exp_maps = sorted(
            (gtc_data_repo / "production_outputs/per_sample_plink_files").glob("*.map")
        )
        assert all([file_hashes_equal(o, e) for o, e in zip(obs_maps, exp_maps)])


def test_create_gtc_merge_list(gtc_data_repo):
    # GIVEN: A RealData repo with the following Snakefile
    gtc_data_repo.make_snakefile(
        """
        from cgr_gwas_qc import load_config

        cfg = load_config()

        include: cfg.rules("entry_points.smk")

        rule all:
            input: "plink_start/mergeList.txt"
        """
    )
    # and the intermediate per_sample_plink_files
    gtc_data_repo.copy("production_outputs/per_sample_plink_files", "per_sample_plink_files")

    # WHEN: I run snakemake
    working_dir = gtc_data_repo.working_dir
    run_snakemake(working_dir, print_stdout=True)

    # THEN: The merge list should contain all Sample_IDs
    n_samples = gtc_data_repo.ss.data.shape[0]
    merge_list = (working_dir / "plink_start/mergeList.txt").read_text().strip().split("\n")
    assert n_samples == len(merge_list)


# def run_and_check_for_plink_start_samples(tmp_path):
#     """Helper function to make sure expected outputs are present."""
#     with chdir(tmp_path):
#         run(["snakemake", "-j1", "--use-conda", "--nocolor"], check=True)
#         assert Path("plink_start/samples.bed").exists()
#         assert Path("plink_start/samples.bim").exists()
#         assert Path("plink_start/samples.fam").exists()


# @pytest.mark.workflow
# def test_gtc_entry(small_gtc_working_dir: Path, conda_envs):
#     """GTC entry point.

#     This is the primary entry point and requires sample level GTC files.
#     """
#     conda_envs.copy_all_envs(small_gtc_working_dir)
#     snake = small_gtc_working_dir / "Snakefile"
#     snake.write_text(snakefile)
#     run_and_check_for_plink_start_samples(small_gtc_working_dir)


# @pytest.mark.workflow
# def test_ped_entry(small_ped_working_dir: Path, conda_envs):
#     """Aggregated PED entry point.

#     This entry point expects an aggregated (multi-sample) PED/MAP file.
#     """
#     conda_envs.copy_all_envs(small_ped_working_dir)
#     snake = small_ped_working_dir / "Snakefile"
#     snake.write_text(snakefile)
#     run_and_check_for_plink_start_samples(small_ped_working_dir)


# @pytest.mark.workflow
# def test_bed_entry(small_bed_working_dir: Path):
#     """Aggregated BED entry point.

#     This entry point expects an aggregated (multi-sample) BED/BIM/FAM file.
#     Since this is the expected output of the entry point module, these files
#     are just symbolically linked to the expected location
#     ``plink_start/samples.{bed,bim,fam}``.
#     """
#     snake = small_bed_working_dir / "Snakefile"
#     snake.write_text(snakefile)
#     run_and_check_for_plink_start_samples(small_bed_working_dir)

#     with chdir(small_bed_working_dir):
#         # Check files are symbolic links
#         assert Path("plink_start/samples.bed").is_symlink()
#         assert Path("plink_start/samples.bim").is_symlink()
#         assert Path("plink_start/samples.fam").is_symlink()
