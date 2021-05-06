import os
import subprocess as sp
import sys
from pathlib import Path
from textwrap import wrap
from typing import List, Optional

import typer
from snakemake.io import load_configfile

from cgr_gwas_qc import load_config
from cgr_gwas_qc.cluster_profiles import env
from cgr_gwas_qc.config import ConfigMgr
from cgr_gwas_qc.version import __version__


def main(
    cgems: bool = typer.Option(False, help="Run using the CGEMs/CCAD cluster profile."),
    biowulf: bool = typer.Option(False, help="Run using the Biowulf cluster profile."),
    cluster_profile: Optional[Path] = typer.Option(None, help="Path to a custom cluster profile."),
    time_hr: int = typer.Option(
        120, help="The walltime limit (in hours) for the main snakemake process."
    ),
    queue: Optional[str] = typer.Option(
        None,
        help="Name of the queue for the main process to use. Only needed if using `cluster_profile`.",
    ),
    submission_cmd: Optional[str] = typer.Option(
        None,
        help="Name of the command to use for submission (i.e., qsub). Only needed if using `cluster_profile`.",
    ),
    dry_run: bool = typer.Option(False, help="Perform a dry-run but don't submit."),
    notemp: bool = typer.Option(False, help="Turn off temporary file deletion."),
):
    """Run CGR GwasQcPipeline on a cluster."""

    check_exclusive_options(cgems, biowulf, cluster_profile)

    payload = {
        "python_executable": sys.executable,
        "working_dir": os.getcwd(),
        "snakefile": ConfigMgr.SNAKEFILE,
        "version": __version__,
        "cgems": cgems,
        "biowulf": biowulf,
        "time_hr": time_hr,
        "local_mem_mb": 1024,
        "local_tasks": 1,
        "added_options": "--notemp " if notemp else "",
    }

    cfg = load_config()
    sample_size = cfg.ss.shape[0]
    if sample_size < 1_000:  # bump up local resources and run a bunch or rules there
        payload["local_tasks"] = 4
        payload["local_mem_mb"] = 1024 * 4
        payload["time_hr"] = 8

    if cgems:
        payload["profile"] = get_profile("cgems")
        payload["queue"] = queue or ("all.q" if time_hr <= 24 else "long.q")
        payload["added_options"] += get_grouping_settings(sample_size)
        submission_cmd = "qsub"
    elif biowulf:
        payload["profile"] = get_profile("biowulf")
        payload["queue"] = queue or ("quick,norm" if time_hr <= 4 else "norm")
        payload["added_options"] += get_grouping_settings(sample_size)
        submission_cmd = "sbatch"
    else:
        payload["profile"] = check_custom_cluster_profile(cluster_profile, queue, submission_cmd)
        payload["queue"] = queue

    run_script = create_submission_script(payload)
    if not dry_run:
        job_id = sp.check_output([submission_cmd, run_script]).decode().strip()  # type: ignore
        print(f"Submitted {job_id}")


def check_exclusive_options(cgems, biowulf, cluster_profile):
    if sum([cgems, biowulf, (cluster_profile is not None)]) > 1:
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
    cgr_profiles = (Path(__file__).parents[1] / "cluster_profiles").resolve()

    if cluster == "cgems":
        return (cgr_profiles / "cgems").as_posix()

    if cluster == "biowulf":
        return (cgr_profiles / "biowulf").as_posix()


def get_group_size(n_samples: int) -> int:
    """Return how many samples to group together"""
    return 500 if n_samples <= 1000 else 1000


def get_per_sample_rules() -> List[str]:
    """Pull list of rules to group from cluster.yml"""
    cluster_file = Path(__file__).parents[1].resolve() / "cluster_profiles/cluster.yaml"
    cluster_profile = load_configfile(cluster_file)
    return sorted([key for key in cluster_profile["group_jobs"].keys()])


def get_grouping_settings(sample_size: int):
    group_size = get_group_size(sample_size)  # 1. Figure out the group size to use
    rules_to_group = get_per_sample_rules()  # 2. Get a list of per sample rules to group

    # 3. Build group setting string
    group_options = []
    if rules_to_group:
        group_options.append("--group-components")
        for rule in rules_to_group:
            group_options.append(f"{rule}={group_size}")

    return " ".join(group_options)


def create_submission_script(payload) -> str:
    template = env.get_template("snakemake.sh")
    run_script = Path(".snakemake/GwasQcPipeline_submission.sh")
    run_script.parent.mkdir(exist_ok=True)
    run_script.write_text(template.render(**payload))
    return run_script.as_posix()
