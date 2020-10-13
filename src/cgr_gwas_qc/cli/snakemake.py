import os
from typing import List, Optional

import typer
from snakemake import snakemake

from cgr_gwas_qc import load_config

app = typer.Typer(add_completion=False)


@app.command(name="snakemake")
def main(
    cores: int = 1,
    nodes: int = 1,
    local_cores: int = 1,
    workdir: str = os.getcwd(),
    dryrun: bool = False,
    touch: bool = False,
    forcetargets: bool = False,
    forceall: bool = False,
    forcerun: List[str] = [],
    printdag: bool = False,
    printrulegraph: bool = False,
    printfilegraph: bool = False,
    printd3dag: bool = False,
    nocolor: bool = False,
    quiet: bool = False,
    keepgoing: bool = False,
    jobname: str = "snakejob.{rulename}.{jobid}.sh",
    lock: bool = True,
    unlock: bool = False,
    conda_cleanup_envs: bool = False,
    debug: bool = False,
    use_conda: bool = False,
    conda_prefix: Optional[str] = None,
    conda_create_envs_only: bool = False,
):
    """Run snakemake on a given snakefile.

    This function provides access to the whole snakemake functionality. It is not thread-safe.

    Args:
        cores (int):                the number of provided cores (ignored when using cluster support) (default 1)
        nodes (int):                the number of provided cluster nodes (ignored without cluster support) (default 1)
        local_cores (int):          the number of provided local cores if in cluster mode (ignored without cluster support) (default 1)
        workdir (str):              path to working directory (default None)
        dryrun (bool):              only dry-run the workflow (default False)
        touch (bool):               only touch all output files if present (default False)
        forcetargets (bool):        force given targets to be re-created (default False)
        forceall (bool):            force all output files to be re-created (default False)
        forcerun (list):            list of files and rules that shall be re-created/re-executed (default [])
        printdag (bool):            print the dag in the graphviz dot language (default False)
        printrulegraph (bool):      print the graph of rules in the graphviz dot language (default False)
        printfilegraph (bool):      print the graph of rules with their input and output files in the graphviz dot language (default False)
        printd3dag (bool):          print a D3.js compatible JSON representation of the DAG (default False)
        nocolor (bool):             do not print colored output (default False)
        quiet (bool):               do not print any default job information (default False)
        keepgoing (bool):           keep goind upon errors (default False)
        jobname (str):              naming scheme for cluster job scripts (default "snakejob.{rulename}.{jobid}.sh")
        lock (bool):                lock the working directory when executing the workflow (default True)
        unlock (bool):              just unlock the working directory (default False)
        conda_cleanup_envs (bool):  just cleanup unused conda environments (default False)
        debug (bool):               allow to use the debugger within rules
        use_conda (bool):           use conda environments for each job (defined with conda directive of rules)
        conda_prefix (str):         the directory in which conda environments will be created (default None)
        conda_create_envs_only (bool):    if specified, only builds the conda environments specified for each job, then exits.
    """

    cfg = load_config()
    snakemake(
        cfg.SNAKEFILE,
        cores=cores,
        nodes=nodes,
        local_cores=local_cores,
        workdir=workdir,
        dryrun=dryrun,
        touch=touch,
        forcetargets=forcetargets,
        forceall=forceall,
        forcerun=forcerun,
        printdag=printdag,
        printrulegraph=printrulegraph,
        printfilegraph=printfilegraph,
        printd3dag=printd3dag,
        nocolor=nocolor,
        quiet=quiet,
        keepgoing=keepgoing,
        jobname=jobname,
        lock=lock,
        unlock=unlock,
        conda_cleanup_envs=conda_cleanup_envs,
        debug=debug,
        use_conda=use_conda,
        conda_prefix=conda_prefix,
        conda_create_envs_only=conda_create_envs_only,
    )


if __name__ == "__main__":
    app()
