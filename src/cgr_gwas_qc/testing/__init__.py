import contextlib
import os
import shutil
from pathlib import Path
from subprocess import PIPE, STDOUT, run
from textwrap import dedent
from typing import Union

from cgr_gwas_qc.config import config_to_yaml
from cgr_gwas_qc.models.config import Config


@contextlib.contextmanager
def chdir(dirname: Union[str, Path]):
    curdir = Path().cwd()
    try:
        os.chdir(dirname)
        yield
    finally:
        os.chdir(curdir)


def make_test_config(current_dir: Path, **kwargs) -> None:
    """Create a GwasQcPipeline config in the working dir"""
    with chdir(current_dir):
        cfg = Config(**kwargs)
        config_to_yaml(cfg)


def make_snakefile(working_dir: Path, contents: str):
    snakefile = working_dir / "Snakefile"
    snakefile.write_text(dedent(contents))


def run_snakemake(working_dir: Path, keep_temp: bool = False, print_stdout: bool = False):
    conda = "mamba" if shutil.which("mamba") else "conda"
    cmd = ["snakemake", "-j1", "--use-conda", "--nocolor", "--conda-frontend", conda]

    if keep_temp:
        cmd.append("--notemp")

    with chdir(working_dir):
        stdout = run(cmd, check=True, stdout=PIPE, stderr=STDOUT).stdout.decode()

    if print_stdout:
        print(stdout)

    return stdout
