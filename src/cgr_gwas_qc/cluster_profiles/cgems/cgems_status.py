#!/usr/bin/env python3

import logging
import re
import shlex
import subprocess as sp
import sys
import time
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
    # if the job has finished it won't appear in qstat and we should check qacct
    # this will also provide the exit status (0 on success, 128 + exit_status on fail)
    # Try getting job with scontrol instead in case sacct is misconfigured
    try:
        qacct_res = sp.check_output(shlex.split(f"qacct -j {job_id}"))
    except sp.CalledProcessError as err:
        logger.warning("qacct process error")
        logger.warning(err)

    match = re.search("exit_status  ([0-9]+)", qacct_res.decode())
    if match:
        exit_code = int(match.group(1))
        return "success" if exit_code == 0 else "failed"

    return None


if __name__ == "__main__":
    main()
