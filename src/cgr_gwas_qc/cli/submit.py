import os
import subprocess as sp
import sys
from pathlib import Path
from textwrap import wrap
from typing import Optional

import typer

from cgr_gwas_qc.cluster_profiles import env


def main(
    cgems: bool = typer.Option(False, help="Run using the CGEMs/CCAD cluster profile."),
    biowulf: bool = typer.Option(False, help="Run using the Biowulf cluster profile."),
    cluster_profile: Optional[Path] = typer.Option(None, help="Path to a custom cluster profile."),
    time_h: int = typer.Option(
        12, help="The walltime limit (in hours) for the main snakemake process."
    ),
    queue: Optional[str] = typer.Option(
        None,
        help="Name of the queue for the main process to use. Only needed if using `cluster_profile`.",
    ),
    submission_cmd: Optional[str] = typer.Option(
        None,
        help="Name of the command to use for submission (i.e., qsub). Only needed if using `cluster_profile`.",
    ),
):
    """Run CGR GwasQcPipeline on a cluster."""

    check_exclusive_options(cgems, biowulf, cluster_profile)

    payload = {
        "python_executable": sys.executable,
        "working_dir": os.getcwd(),
        "cgems": cgems,
        "biowulf": biowulf,
        "time_h": time_h,
    }

    if cgems:
        payload["profile"] = get_profile("cgems")
        payload["queue"] = queue or ("all.q" if time_h <= 24 else "long.q")
        submission_cmd = "qsub"
    elif biowulf:
        payload["profile"] = get_profile("biowulf")
        payload["queue"] = queue or ("all.q" if time_h <= 24 else "long.q")
        submission_cmd = "sbatch"
    else:
        payload["profile"] = check_custom_cluster_profile(cluster_profile, queue, submission_cmd)
        payload["queue"] = queue

    run_script = create_submission_script(payload)
    sp.run([submission_cmd, run_script], shell=True, check=True)  # type: ignore


def check_exclusive_options(cgems, biowulf, cluster_profile):
    if not (cgems ^ biowulf ^ (cluster_profile is not None)):
        typer.echo(
            "\n".join(
                wrap(
                    "Please only provide one of `--cgems`, `--biowulf`, or `--cluster_profile`. "
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
    cgr_profiles = (Path(__file__).parents[2] / "cluster_profiles").resolve()

    if cluster == "cgems":
        return (cgr_profiles / "cgems").as_posix()

    if cluster == "biowulf":
        return (cgr_profiles / "biowulf").as_posix()


def create_submission_script(payload) -> Path:
    template = env.get_template("snakemake.sh")
    run_script = Path(".snakemake/GwasQcPipeline_submission.sh")
    run_script.parent.mkdir(exist_ok=True)
    run_script.write_text(template.render(**payload))
    return run_script
