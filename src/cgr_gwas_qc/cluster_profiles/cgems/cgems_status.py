#!/usr/bin/env python3
import logging
import shlex
import subprocess as sp
import sys
import time
from pathlib import Path
from typing import Optional

logger = logging.getLogger("__name__")
logger.setLevel(40)

MAX_STATUS_ATTEMPTS = 20


def main():
    job_id = int(sys.argv[1])
    for _ in range(MAX_STATUS_ATTEMPTS):
        job_status = check_queue(job_id) or check_job_history(job_id)
        if job_status:
            break
        time.sleep(5)

    print(job_status or "failed")


def check_queue(job_id: int) -> Optional[str]:
    try:
        qstat_res = sp.check_output(shlex.split("qstat -s pr")).decode().strip()
    except sp.CalledProcessError as err:
        logger.error("qstat process error")
        logger.error(err)

    try:
        queue_status = {int(x.split()[0]): x.split()[4] for x in qstat_res.splitlines()[2:]}
        return "failed" if "E" in queue_status[job_id] else "running"
    except KeyError:
        logger.info(f"{job_id} no in queue")

    return None


def check_job_history(job_id: int) -> Optional[str]:
    # If the job has finished it won't appear in qstat. Typically, I would use
    # qacct to check the exit code for each job_id, but it is extremely slow.
    # Instead I am parsing each job's log and making sure it is marked as
    # (100%) done.
    try:
        log = next(Path("logs").glob(f"*.{job_id}")).read_text().strip()
        return "success" if log.endswith("(100%) done") else None
    except StopIteration:
        return None


if __name__ == "__main__":
    main()
