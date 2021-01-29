import os

import typer
from snakemake import snakemake

from cgr_gwas_qc import load_config
from cgr_gwas_qc.config import ConfigMgr

app = typer.Typer(add_completion=False)


def fmt_doc_string(func):
    """Hack to add full path to snakefile in the doc string."""
    func.__doc__ = func.__doc__.format(ConfigMgr.SNAKEFILE)
    return func


@app.command(name="snakemake")
@fmt_doc_string
def main(
    cores: int = typer.Option(
        1,
        "--cores",
        "--jobs",
        "-j",
        help="the number of provided cores (ignored when using cluster support)",
    ),
    local_cores: int = typer.Option(
        1,
        "--local-cores",
        help="the number of provided local cores if in cluster mode (ignored without cluster support)",
    ),
    dryrun: bool = typer.Option(
        False, "--dry-run", "--dryrun", "-n", help="only dry-run the workflow"
    ),
    touch: bool = typer.Option(
        False, "--touch", "-t", help="only touch all output files if present"
    ),
    forcetargets: bool = typer.Option(
        False, "--force", "-f", help="force given targets to be re-created"
    ),
    forceall: bool = typer.Option(
        False, "--forceall", "-F", help="force all output files to be re-created"
    ),
    quiet: bool = typer.Option(
        False, "--quiet", "-q", help="do not print any default job information"
    ),
    keepgoing: bool = typer.Option(False, "--keep-going", "-k", help=" keep going upon errors"),
    unlock: bool = typer.Option(False, help=" just unlock the working directory"),
    conda_cleanup_envs: bool = typer.Option(
        False, "--conda-cleanup-envs", help="just cleanup unused conda environments"
    ),
    use_conda: bool = typer.Option(
        False, "--use-conda", help="use conda environments for each job"
    ),
    conda_create_envs_only: bool = typer.Option(
        False,
        "--conda-create-envs-only",
        help="if specified, only builds the conda environments specified for each job, then exits.",
    ),
):
    """Run the Gwas Qc Pipeline using Snakemake.

    We only expose a few common snakemake options. If you want to run
    snakemake directly then use:

        snakemake -s {} OPTIONS TARGETS
   """

    cfg = load_config()
    workdir = os.getcwd()

    snakemake(
        cfg.SNAKEFILE,
        cores=cores,
        local_cores=local_cores,
        workdir=workdir,
        dryrun=dryrun,
        touch=touch,
        forcetargets=forcetargets,
        forceall=forceall,
        quiet=quiet,
        keepgoing=keepgoing,
        unlock=unlock,
        conda_cleanup_envs=conda_cleanup_envs,
        use_conda=use_conda,
        conda_create_envs_only=conda_create_envs_only,
    )


if __name__ == "__main__":
    app()
