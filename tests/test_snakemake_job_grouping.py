import shutil
import subprocess as sp
from textwrap import dedent

import pytest
from snakemake.utils import read_job_properties

from cgr_gwas_qc.testing import chdir, make_snakefile
from cgr_gwas_qc.testing.data import FakeData


def test_basic_grouping(tmp_path, qsub):
    """Check job properties of grouped jobs.

    Here is the most basic job group. The job script just runs for the output
    of `b/1.out`.
    """
    shutil.copyfile(FakeData() / "job_scripts/basic_group.smk", tmp_path / "Snakefile")
    with chdir(tmp_path):
        sp.run(["snakemake", "--cluster", "qsub", "-j", "1"], check=True)

    expected_properties = {
        "cluster": {},
        "groupid": "grp0",
        "input": [],
        "local": False,
        "output": ["b/1.out"],
        "resources": {},
        "threads": 1,
        "type": "group",
    }

    job_properties = read_job_properties(tmp_path / "job_script.sh")
    del job_properties["jobid"]  # this is random so ignore

    assert expected_properties == job_properties


def test_wildcard_grouping(tmp_path, qsub):
    """Curious if you can use wildcards to group rules.

    Yes, it looks like if you use a wildcard it will just replace that as the
    group name. This will allow general rules, that operate in multiple
    places, to be grouped with their upstream components. I am thinking about
    `plink --maf` which runs in multiple places, so it would be nice to group
    it since it is a fast step.
    """
    make_snakefile(
        tmp_path,
        dedent(
            """
            rule all:
                input: expand("b/{sample}.out", sample=[1])

            rule a:
                output: "a/{sample}.out"
                group: "{sample}"
                shell: "touch {output}"

            rule b:
                input: "a/{sample}.out"
                output: "b/{sample}.out"
                group: "{sample}"
                shell: "touch {output}"
            """
        ),
    )
    with chdir(tmp_path):
        sp.run(["snakemake", "--cluster", "qsub", "-j", "1"], check=True)

    expected_properties = {
        "cluster": {},
        "groupid": "1",  # group-id is set as the value of {sample}
        "input": [],
        "local": False,
        "output": ["b/1.out"],
        "resources": {},
        "threads": 1,
        "type": "group",
    }

    job_properties = read_job_properties(tmp_path / "job_script.sh")
    del job_properties["jobid"]  # this is random so ignore

    assert expected_properties == job_properties


@pytest.mark.parametrize("a_mem,b_mem", [(2, 1), (1, 2)])
def test_basic_grouping_max_resource(a_mem, b_mem, tmp_path, qsub):
    """How do groups work with resources?

    Groups appear to take the maximum resource value for all resources. So if
    a resource is defined in two rules it will use the max. If it is defined
    in only one rule it will still use that resource.
    """
    make_snakefile(
        tmp_path,
        dedent(
            f"""
            rule all:
                input: expand("b/{{sample}}.out", sample=[1])

            rule a:
                output: "a/{{sample}}.out"
                group: 0
                resources:
                    mem={a_mem},
                    ssd=2
                shell: "touch {{output}}"

            rule b:
                input: "a/{{sample}}.out"
                output: "b/{{sample}}.out"
                group: 0
                resources:
                    mem={b_mem},
                shell: "touch {{output}}"
            """
        ),
    )
    with chdir(tmp_path):
        sp.run(["snakemake", "--cluster", "qsub", "-j", "1"], check=True)

    expected_properties = {
        "cluster": {},
        "groupid": 0,
        "input": [],
        "local": False,
        "output": ["b/1.out"],
        "resources": {"mem": 2, "ssd": 2},  # max mem is always used
        "threads": 1,
        "type": "group",
    }

    job_properties = read_job_properties(tmp_path / "job_script.sh")
    del job_properties["jobid"]  # this is random so ignore

    assert expected_properties == job_properties


@pytest.mark.parametrize("a_cpu,b_cpu", [(2, 1), (1, 2)])
def test_basic_grouping_max_threads(a_cpu, b_cpu, tmp_path, qsub):
    """How do groups work with threads?

    Theads behave the same as resources, the max value is used from all rules
    in the group.
    """
    make_snakefile(
        tmp_path,
        dedent(
            f"""
            rule all:
                input: expand("b/{{sample}}.out", sample=[1])

            rule a:
                output: "a/{{sample}}.out"
                group: 0
                threads: {a_cpu}
                shell: "touch {{output}}"

            rule b:
                input: "a/{{sample}}.out"
                output: "b/{{sample}}.out"
                group: 0
                threads: {b_cpu}
                shell: "touch {{output}}"
            """
        ),
    )
    with chdir(tmp_path):
        sp.run(["snakemake", "--cluster", "qsub", "-j", "1"], check=True)

    expected_properties = {
        "cluster": {},
        "groupid": 0,
        "input": [],
        "local": False,
        "output": ["b/1.out"],
        "resources": {},
        "threads": 2,  # should be the max
        "type": "group",
    }

    job_properties = read_job_properties(tmp_path / "job_script.sh")
    del job_properties["jobid"]  # this is random so ignore

    assert expected_properties == job_properties


def test_basic_component_grouping(tmp_path, qsub):
    """How does component grouping work?

    Snakemake as `--group-components` which can group multiple runs of the
    same rule together.

    .. warning::
        Component grouping sums the number of threads and resources.
    """
    shutil.copyfile(
        FakeData() / "job_scripts/single_rule_10_samples_with_resources.smk", tmp_path / "Snakefile"
    )
    with chdir(tmp_path):
        sp.run(
            ["snakemake", "--cluster", "qsub", "-j", "1", "--group-components", "grp0=5"],
            check=True,
        )

    job_properties = read_job_properties(tmp_path / "job_script.sh")

    assert len(job_properties["output"]) == 5  # 5 samples per job
    assert job_properties["threads"] == 10  # NOTE: threads is the sum!!
    assert job_properties["resources"]["mem"] == 10  # NOTE: resources use the sum!!


def test_basic_component_grouping_cluster_override(tmp_path, qsub):
    """How does component grouping work?

    Snakemake as `--group-components` which can group multiple runs of the
    same rule together.

    .. warning::
        Component grouping sums the number of threads and resources.
    """
    make_snakefile(
        tmp_path,
        dedent(
            """
            rule all:
                input: expand("b/{sample}.out", sample=[x for x in range(1, 11)])

            rule b:
                output: "b/{sample}.out"
                group: "grp0"
                threads: 1
                resources:
                    mem=2
                shell: "touch {output}"
            """
        ),
    )
    with chdir(tmp_path):
        sp.run(
            ["snakemake", "--cluster", "qsub", "-j", "2", "--group-components", "grp0=5"],
            check=True,
        )

    job_properties = read_job_properties(tmp_path / "job_script.sh")

    assert len(job_properties["output"]) == 5  # 5 samples per job
    assert job_properties["threads"] == 5  # NOTE: threads is the sum!!
