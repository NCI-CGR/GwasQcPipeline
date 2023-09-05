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

STATUS_CODES = {
    "BOOT_FAIL": "failed",
    "CANCELLED": "failed",
    "COMPLETED": "success",
    "DEADLINE": "failed",
    "FAILED": "failed",
    "NODE_FAIL": "failed",
    "OUT_OF_MEMORY": "failed",
    "PENDING": "running",
    "PREEMPTED": "failed",
    "RUNNING": "running",
    "REQUEUED": "running",
    "RESIZING": "running",
    "REVOKED": "running",
    "SUSPENDED": "failed",
    "TIMEOUT": "failed",
}

def main():
    job_id = int(sys.argv[1])
    for _ in range(MAX_STATUS_ATTEMPTS):
        job_status = check_sacct(job_id) or check_scontrol(job_id)
        if job_status:
            break
        time.sleep(5)

    print(job_status or "failed")


def check_sacct(job_id: int) -> Optional[str]:
    try:
        job_info = sp.check_output(shlex.split(f"sacct -P -b -j {job_id} -n"))      
    except sp.CalledProcessError as err:
        logger.error("sacct process error")
        logger.error(err)
        return None

    try:
        status = {x.split("|")[0]: x.split("|")[1] for x in job_info.decode().strip().split("\n")}
        return STATUS_CODES.get(status[f"{job_id}"], None)
    except IndexError:
        return None


def check_scontrol(job_id: int) -> Optional[str]:
    try:
        job_info = sp.check_output(shlex.split(f"scontrol -o show job {job_id}"))
    except sp.CalledProcessError as err:
        logger.error("scontrol process error")
        logger.error(err)
        return None

    m = re.search(r"JobState=(\w+)", job_info.decode())
    status = {job_id: m.group(1)} if m else {}
    return STATUS_CODES.get(status[job_id], None)


if __name__ == "__main__":
    main()
