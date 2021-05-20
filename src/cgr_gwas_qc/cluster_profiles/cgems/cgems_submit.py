#!/usr/bin/env python3
from dataclasses import dataclass, field
from typing import List

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
class CgemsOptions(ClusterOptions):
    queue: List[str] = field(
        default_factory=lambda: [
            "all.q",
            "seq-alignment.q",
            "seq-calling.q",
            "seq-calling2.q",
            "seq-gvcf.q",
        ]
    )
    log: str = "logs/{rulename}_{job_id}.$JOB_ID"

    def __str__(self):
        # See cgems_jobscript.sh for default sge options
        cmd = "qsub -q {queue} -N gqc.{job_id} -o {log}"

        if self.threads > 1:
            cmd += " -pe by_node {threads}"

        return cmd.format(
            queue=",".join(self.queue),
            mem=self.mem_mb,
            time=str(self.time),
            threads=self.threads,
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

    cgems_options = CgemsOptions.from_dict(options)

    # Submit job script
    cmd = cgems_options.command + [jobscript]
    print(submit_job(cmd))


if __name__ == "__main__":
    main()
