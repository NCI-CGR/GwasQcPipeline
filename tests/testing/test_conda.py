from cgr_gwas_qc.conda import get_snakemake_conda_env
from cgr_gwas_qc.testing import make_snakefile, run_snakemake
from cgr_gwas_qc.testing.conda import CondaEnv


def test_copy_conda_env(tmp_path):
    """Test copying conda environment to another folder.

    This should build a cached version of each conda environment and then
    copy them into the working directory for use with snakemake.
    """
    # GIVEN: our conda object
    conda_env = CondaEnv()
    # and we get the correct renamed version of the environment
    env_file_name = conda_env.conda_env_path / "plink2.yml"
    local_copy = get_snakemake_conda_env(env_file_name, tmp_path)

    # WHEN: we copy an environment into a temporary directory
    conda_env.copy_env("plink2", tmp_path)

    # THEN: the renamed environment
    # should exist
    assert local_copy.exists()
    # and have the plink binary
    assert (local_copy / "bin/plink").exists()


def test_conda_env_works_with_snakemake(tmp_path):
    """Will conda env copies work with Snakemake"""
    # GIVEN: a conda environment being called in a snakemake
    conda_env = CondaEnv()
    env_file_name = conda_env.conda_env_path / "plink2.yml"
    make_snakefile(
        tmp_path,
        f"""
        rule plink_test:
            output: "test.txt"
            conda: "{env_file_name}"
            shell: "plink --help > {{output[0]}}"
        """,
    )

    # WHEN: we copy the environment to the working directory and run snakemake
    conda_env.copy_env("plink2", tmp_path)
    stdout = run_snakemake(tmp_path)

    # THEN:
    # snakemake should not have to create conda environment
    assert "Creating conda environment" not in stdout
    # the output should be created
    assert (tmp_path / "test.txt").exists()
    # the plink help message should contain the word plink
    assert "plink" in (tmp_path / "test.txt").read_text().lower()
