#!/usr/bin/env python3
from dataclasses import dataclass, field
from datetime import timedelta
from typing import Set

from snakemake.utils import read_job_properties

from cgr_gwas_qc.cluster_profiles import (
    ClusterOptions,
    load_cluster_config,
    parse_jobscript,
    submit_job,
    update_group_properties,
    update_properties,
)


@dataclass
class BiowulfOptions(ClusterOptions):
    queue: Set[str] = field(default_factory=lambda: {"quick", "norm"})
    log: str = "logs/{rulename}_{job_id}.%j"

    def __str__(self):
        # See cgems_jobscript.sh for default sge options
        cmd = "sbatch --partition={queue} --job-name={rulename}.{job_id} --cpus-per-task={threads} --mem={mem_gb:0.0f}gb --time={time} --output={log}"

        if self.time > timedelta(hours=4):
            self.queue.discard("quick")

        return cmd.format(
            queue=",".join(self.queue),
            mem_gb=self.mem_gb,
            time=str(self.time),
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

    # Build cgems options
    update_properties(options, cluster_config, job_properties)

    if job_properties["type"] == "group":
        update_group_properties(options, cluster_config, job_properties, jobscript)

    biowulf_options = BiowulfOptions.from_dict(options)

    # Submit job script
    cmd = biowulf_options.command + [jobscript]
    print(submit_job(cmd))


if __name__ == "__main__":
    main()
