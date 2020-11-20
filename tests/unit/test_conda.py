from pathlib import Path
from subprocess import run
from textwrap import dedent

from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.testing.conda import CondaEnv


def test_copy_conda_env(tmp_path):
    """Test copying conda environment to another folder.

    This should build a cached version of each conda environment and then
    copy them into the working directory for use with snakemake.
    """
    conda_env = CondaEnv()
    conda_env.copy_env("plink2", tmp_path)

    env_file_name = conda_env.conda_env_path / "plink2.yml"
    local_copy = conda_env._get_working_dir_env_path(env_file_name, tmp_path)

    assert local_copy.exists()
    assert (local_copy / "bin/plink").exists()


def test_conda_env_works_with_snakemake(tmp_path):
    """Will conda env copies work with Snakemake"""
    conda_env = CondaEnv()
    conda_env.copy_env("plink2", tmp_path)

    env_file_name = conda_env.conda_env_path / "plink2.yml"

    snakefile = dedent(
        f"""\
        rule plink_test:
            output: "test.txt"
            conda: "{env_file_name}"
            shell: "plink --help > {{output[0]}}"
        """
    )

    snake = tmp_path / "Snakefile"
    snake.write_text(snakefile)

    with chdir(tmp_path):
        run(["snakemake", "-j1", "--use-conda", "--nocolor"], check=True)
        assert Path("test.txt").exists()
        assert "plink" in Path("test.txt").read_text().lower()
