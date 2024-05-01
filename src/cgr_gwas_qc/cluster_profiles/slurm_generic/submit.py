#!/usr/bin/env python3

from dataclasses import dataclass  # , field

from snakemake.utils import read_job_properties

from cgr_gwas_qc import load_config
from cgr_gwas_qc.cluster_profiles import (
    ClusterOptions,
    load_cluster_config,
    parse_jobscript,
    submit_job,
    update_group_properties,
    update_properties,
)

# from datetime import timedelta
# from typing import Set


cfg = load_config()
queue = cfg.config.slurm_partition


@dataclass
class SlurmOptions(ClusterOptions):
    log: str = "logs/{rulename}_{job_id}.%j"

    def __str__(self):
        cmd = (
            "sbatch"
            " --partition={queue}"
            " --job-name={rulename}.{job_id}"
            " --cpus-per-task={threads}"
            " --mem={mem_gb:0.0f}gb"
            " --time={time}"
            " --output={log}"
            " --parsable"
        )

        formatted_time = f"{self.time.days}-{self.time.seconds // 3600}"  # days-hours for slurm

        return cmd.format(
            queue=queue,
            mem_gb=self.mem_gb,
            time=formatted_time,
            threads=self.threads,
            rulename=self.rulename,
            job_id=self.job_id,
            log=self.log.format(rulename=self.rulename, job_id=self.job_id),
        )


def main(jobscript=None):
    # Load default options from cluster.yml
    cluster_config = load_cluster_config()
    options = cluster_config.get("__default__", {}).copy()

    # Load the job script
    jobscript = jobscript or parse_jobscript()
    job_properties = read_job_properties(jobscript)

    # Build slurm options
    update_properties(options, cluster_config, job_properties)

    if job_properties["type"] == "group":
        update_group_properties(options, cluster_config, job_properties, jobscript)

    slurm_options = SlurmOptions.from_dict(options)

    # Submit job script
    cmd = slurm_options.command + [jobscript]
    print(submit_job(cmd))


if __name__ == "__main__":
    main()
