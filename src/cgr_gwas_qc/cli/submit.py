import os
import subprocess as sp
import sys
from pathlib import Path
from textwrap import wrap
from typing import Optional

import typer

from cgr_gwas_qc import load_config
from cgr_gwas_qc.cluster_profiles import env
from cgr_gwas_qc.config import ConfigMgr
from cgr_gwas_qc.version import __version__

app = typer.Typer(add_completion=False)


@app.command()
def main(
    cgems: bool = typer.Option(
        False, help="Run the workflow using the CGEMs/CCAD cluster profile."
    ),
    biowulf: bool = typer.Option(
        False, help="Run the workflow using the cli Biowulf cluster profile."
    ),
    ccad2: bool = typer.Option(False, help="Run the workflow using the cli ccad2 cluster profile."),
    slurm: bool = typer.Option(
        False,
        help="Run the workflow using the generic slurm cluster profile. "
        "This option requires that you specify the ``slurm_partition`` "
        "in the ``config.yml`` file.",
    ),
    cluster_profile: Optional[Path] = typer.Option(
        None,
        help="Path to a custom cluster profile. See https://github.com/snakemake-profiles/doc.",
    ),
    subworkflow: Optional[str] = typer.Option(
        None, help="Run a sub-workflow instead of the full workflow."
    ),
    time_hr: int = typer.Option(
        120,
        help="The walltime limit (in hours) for the main snakemake process. "
        "Ignored if using ``--cgems``.",
    ),
    queue: Optional[str] = typer.Option(
        None,
        help="Name of the queue for the main process to use. "
        "This option is required if you are using ``--cluster-profile``. "
        "CGEMs/CCAD, CCAD2, and Biowulf users may want to use this to override which queue to use for the main process.",
    ),
    submission_cmd: Optional[str] = typer.Option(
        None,
        help="Name of the command to use for submission (i.e., qsub, sbatch). "
        "This is required if you are using ``--cluster-profile``. "
        "Ignored if using ``--cgems`` or ``--biowulf`` or ``--ccad2``.",
    ),
    dry_run: bool = typer.Option(
        False,
        help="Create the submission script but do not send to the schedular. "
        "The generate submission script is in ``.snakemake/GwasQcPipeline_submission.sh``. "
        "This file can be edited and submitted directly to the cluster.",
    ),
    notemp: bool = typer.Option(
        False,
        help="Do not delete temporary files. "
        "This option runs the workflow using ``snakemake --notemp``.",
    ),
    local_mem_mb: int = typer.Option(
        8000,
        help="The amount of memory to use for main snakemake process and local rules [default 8GB]. "
        "There are a number of rules that are run along side of the main snakemake process instead of submitting a new job. "
        "If you have a very large project >50k samples, you may want to increase the amount of memory used for this job. "
        "Ignored if using ``--cgems``.",
    ),
    local_tasks: int = typer.Option(
        4,
        help="The number of threads for the main snakemake process and local rules."
        "There are a number of rules that are run along side of the main snakemake process instead of submitting a new job. "
        "If you have a very large project >50k samples, you may want to increase the number of CPUs used for this job. "
        "Ignored if using ``--cgems``.",
    ),
):
    """Submit the CGR GwasQcPipeline to a cluster for execution.

    The ``cgr submit`` command will create a submission script and submit it to the cluster.
    We will create an optimized submission script for users of CGEMs/CCAD, CCAD2 and Biowulf.
    For other systems you will need to provide a snakemake cluster profile to tell snakemake how use your system.
    If you are submitting to a Slurm cluster, you can use the generic Slurm cluster profile we provide by including the ``--slurm`` option in your submit command.

    Users running on CGEMs/CCAD will typically run::

        cgr submit --cgems

    Users running on Biowulf will typically run::

        cgr submit --biowulf

    Users running on CCAD2 will typically run::

        cgr submit --ccad2

    Users running on other Slurm systems will typically run::

        cgr submit --slurm

    Users running with a custom cluster profile will typically run:

    .. code-block:: none

        cgr submit \\
            --profile <path to custom cluster profile> \\
            --queue <name of the queue to submit main job> \\
            --submission-cmd <tool used for submit such as sbatch or qsub>

    .. note::
        Sometimes it may be useful to edit the submission script before submitting it to the cluster.
        In that case you can add the ``--dry-run`` option and edit the generated file in ``.snakemake/GwasQcPipeline_submission.sh``.
        You would then submit this script directly to your cluster (i.e., ``qsub .snakemake/GwasQcPipeline_submission.sh``).
    """

    check_exclusive_options(cgems, biowulf, ccad2, slurm, cluster_profile)

    payload = {
        "python_executable": sys.executable,
        "working_dir": os.getcwd(),
        "snakefile": ConfigMgr.SNAKEFILE,
        "version": __version__,
        "cgems": cgems,
        "biowulf": biowulf,
        "ccad2": ccad2,
        "slurm": slurm,
        "time_hr": time_hr,
        "local_mem_mb": local_mem_mb,
        "local_tasks": local_tasks,
        "added_options": "",
    }
    snake_config = {"cluster_mode": True, "notemp": False}

    if notemp:
        payload["added_options"] += "--notemp "  # type: ignore
        snake_config["notemp"] = True

    if subworkflow:
        payload["added_options"] += f"--subworkflow {subworkflow} "  # type: ignore

    # add global config options only used in cluster mode.
    payload["added_options"] += "--config {} ".format(  # type: ignore
        " ".join("=".join([k, str(v)]) for k, v in snake_config.items())
    )

    cfg = load_config()
    sample_size = cfg.ss.shape[0]
    if sample_size < 1_000:  # need less walltime for smaller sample size
        payload["time_hr"] = 8

    if cgems:
        payload["profile"] = get_profile("cgems")
        payload["queue"] = queue or ("all.q" if time_hr <= 24 else "long.q")
        submission_cmd = "qsub"
    elif biowulf:
        payload["profile"] = get_profile("biowulf")
        payload["queue"] = queue or ("quick,norm" if time_hr <= 4 else "norm")
        submission_cmd = "sbatch"
    elif ccad2:
        payload["profile"] = get_profile("ccad2")
        payload["queue"] = queue or ("defq,defq" if time_hr <= 4 else "defq")
        submission_cmd = "sbatch"
    elif slurm:
        payload["profile"] = get_profile("slurm")
        queue = cfg.config.slurm_partition
        if queue is None:
            raise ValueError(
                """
            \033[1;31mMissing required configuration key-value pair for slurm option\033[0m\n
            The `--slurm` option requires the `slurm_partition`
            key-value pair to be defined in your `config.yml` file.
            This pair specifies the partition where your job will
            run on the Slurm scheduler. The key is `slurm_partition`
            and the value should be the name of the desired partition.
            You can execute `sinfo` to see which partitions are available.
            Add a line like this in your config.yml:\n
                    \033[32mslurm_partition: <partition_name>\033[0m
            """
            )
        payload["queue"] = queue
        submission_cmd = "sbatch"
    else:
        payload["profile"] = check_custom_cluster_profile(cluster_profile, queue, submission_cmd)
        payload["queue"] = queue

    run_script = create_submission_script(payload)

    if not dry_run:
        job_id = sp.check_output([submission_cmd, run_script]).decode().strip()  # type: ignore
        print(f"Submitted {job_id}")


def check_exclusive_options(cgems, biowulf, ccad2, slurm, cluster_profile):
    if sum([cgems, biowulf, ccad2, slurm, (cluster_profile is not None)]) > 1:
        typer.echo(
            "\n".join(
                wrap(
                    "Please only provide one of `--cgems`, `--biowulf`, `--ccad2, `--slurm`, or `--cluster_profile`. "
                    "Run `cgr submit --help` for more information.",
                    width=100,
                )
            )
        )
        raise typer.Exit(code=1)


def check_custom_cluster_profile(cluster_profile, queue, submission_cmd):
    if queue is None:
        raise ValueError("You must provide a queue to use with `cluster_profile`.")

    if submission_cmd is None:
        raise ValueError("You must provide a submission_cmd with `cluster_profile`.")

    if not cluster_profile.exists() or not cluster_profile.is_dir():
        raise ValueError(
            "\n".join(
                wrap(
                    "`cluster_profile` should point to a directory with a snakemake "
                    "cluster profile as described at https://github.com/Snakemake-Profiles/doc",
                    width=100,
                )
            )
        )

    return cluster_profile.resolve().as_posix()


def get_profile(cluster: str):
    cgr_profiles = (Path(__file__).parents[1] / "cluster_profiles").resolve()

    if cluster == "cgems":
        return (cgr_profiles / "cgems").as_posix()

    if cluster == "biowulf":
        return (cgr_profiles / "biowulf").as_posix()

    if cluster == "ccad2":
        return (cgr_profiles / "ccad2").as_posix()

    if cluster == "slurm":
        return (cgr_profiles / "slurm_generic").as_posix()
        # if queue is None:
        #    raise ValueError("You must provide provide slurm_parition in the config.yml to use with `--slurm`.")


def create_submission_script(payload) -> str:
    template = env.get_template("snakemake.sh")
    run_script = Path(".snakemake/GwasQcPipeline_submission.sh")
    run_script.parent.mkdir(exist_ok=True)
    run_script.write_text(template.render(**payload))
    return run_script.as_posix()


if __name__ == "__main__":
    app()

typer_click_object = typer.main.get_command(app)  # only needed for building documentation
