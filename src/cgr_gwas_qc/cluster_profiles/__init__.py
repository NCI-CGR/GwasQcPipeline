import argparse
import logging
import re
import subprocess as sp
from dataclasses import dataclass
from datetime import timedelta
from pathlib import Path
from typing import Dict, List

from jinja2 import Environment, PackageLoader
from snakemake import io

logger = logging.getLogger("__name__")
logger.setLevel(40)


env = Environment(
    loader=PackageLoader("cgr_gwas_qc", "cluster_profiles"),
    trim_blocks=True,
    lstrip_blocks=True,
    keep_trailing_newline=True,
)


@dataclass
class ClusterOptions:
    rulename: str
    job_id: int
    threads: int
    mem: int
    time: timedelta

    @classmethod
    def from_dict(cls, data: Dict):
        return cls(
            data["rulename"],
            data["job_id"],
            data["threads"],
            cls._get_mem(data),
            cls._get_time(data),
        )

    @property
    def mem_kb(self):
        return self.mem / 1024

    @property
    def mem_mb(self):
        return self.mem / 1024 ** 2

    @property
    def mem_gb(self):
        return self.mem / 1024 ** 3

    @property
    def command(self):
        return str(self).split()

    @staticmethod
    def _get_mem(data: Dict):
        if "mem_kb" in data:
            return 1024 * int(data["mem_kb"])  # convert KB to Bytes

        if "mem_mb" in data:
            return 1024 ** 2 * int(data["mem_mb"])  # convert MB to Bytes

        if "mem_gb" in data:
            return 1024 ** 3 * int(data["mem_gb"])  # Convert GB to Bytes

        if "mem" in data:  # Assume mem is mem_gb
            return 1024 ** 3 * int(data["mem"])  # Convert GB to Bytes

        return 1024 ** 3 * 2  # Default to 2GB

    @staticmethod
    def _get_time(data: Dict):
        if "time_min" in data:
            return timedelta(minutes=data["time_min"])

        if "time_hr" in data:
            return timedelta(hours=data["time_hr"])

        if "time" in data:  # Assume time is time_hr
            return timedelta(hours=data["time"])

        return timedelta(hours=1)  # Default to 1 hr

    def __post_init__(self):
        """Ensure that the log dir exists."""
        if hasattr(self, "log"):
            Path(self.log).parent.mkdir(exist_ok=True, parents=True)  # type: ignore # noqa

    def __str__(self):
        raise NotImplementedError


def load_cluster_config(cluster_file=None) -> Dict:
    """Load config to dict either from absolute path or relative to profile dir."""
    if cluster_file is None:
        cluster_file = Path(__file__).parent.resolve() / "cluster.yaml"
    config = io.load_configfile(cluster_file)

    if "__default__" not in config:
        config["__default__"] = {}

    return config


def parse_jobscript() -> str:
    """Minimal CLI to require/only accept single positional argument."""
    p = argparse.ArgumentParser(description="SGE snakemake submit script")
    p.add_argument("jobscript", help="Snakemake jobscript with job properties.")
    return p.parse_args().jobscript


def update_properties(options: Dict, cluster_config: Dict, job_properties: Dict) -> None:
    # Bull basic properties from job script
    options["rulename"] = job_properties.get("rule")
    options["job_id"] = job_properties.get("jobid")
    options["threads"] = job_properties.get("threads")

    # Load resources
    _update_cluster_options(options, cluster_config.get(options["rulename"], {}))  # cluster.yml
    options.update(job_properties.get("cluster", {}))  # --cluster-config
    options.update(job_properties.get("resources", {}))  # resources directive

    return None


def _remove_time(options: Dict):
    keys = list(options.keys())
    for key in keys:
        if key.startswith("time"):
            del options[key]


def _remove_mem(options: Dict):
    keys = list(options.keys())
    for key in keys:
        if key.startswith("mem"):
            del options[key]


def _update_cluster_options(options: Dict, new_options: Dict):
    for key, value in new_options.items():
        if key.startswith("time"):
            _remove_time(options)

        if key.startswith("mem"):
            _remove_mem(options)

        options[key] = value


def update_group_properties(
    options: Dict, cluster_config: Dict, job_properties: Dict, jobscript: str
):
    update_properties(options, cluster_config, job_properties)
    rulenames = _get_rule_names(jobscript)
    n_rules = len(rulenames)

    unique_rulenames = set(rulenames)
    n_unique_rules = len(unique_rulenames)

    rulename = ".".join(sorted(unique_rulenames))  # concatenate rulenames: rule1.rule2
    options["rulename"] = f"GROUP.{rulename[:100]}".rstrip(".")

    # Adjust multi-sample group resources to sane values. NOTE: this only will
    # work well if the cluster correctly use resource masks to limit the number
    # of CPUs.
    group_jobs = cluster_config.get("group_jobs", {})
    if (n_unique_rules == 1) & ((_rulename := rulenames[0]) in group_jobs):
        n_samples = n_rules

        if n_samples < 1000:
            rule_options = group_jobs[_rulename]["small"]
        else:
            rule_options = group_jobs[_rulename]["large"]

        _update_cluster_options(options, rule_options)

    return None


def submit_job(cmd: List[str]):
    return sp.check_output(cmd).decode().rstrip()


def _get_rule_names(jobscript: str) -> List[str]:
    _jobscript = Path(jobscript).read_text()
    match = re.search(r"--allowed-rules (.*?) --", _jobscript)

    if match is None:
        raise ValueError(f"Grouped jobscript missing `--allowed-rules`.\n{_jobscript}")

    allowed_rules = match.group(1)
    return allowed_rules.strip().split()
