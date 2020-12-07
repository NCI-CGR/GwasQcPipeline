import re
from pathlib import Path
from typing import Tuple

import pytest

from cgr_gwas_qc.testing import file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import RealData


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_split_call_rate_filters(tmp_path, conda_envs):
    """Check call rate filter order.

    The legacy workflow runs marker and sample call rate filters
    simultaneously. According to the logs it looks like plink actually runs
    samples then marker filters sequentially. We decided to split this into
    two steps and reverse the order (i.e., snp then sample filtering). Here I
    am testing the splitting part to make sure we get the same results. In
    other words, I am keep in the legacy workflow order (sample -> snp) and
    to make sure splitting does not change anything.
    """
    # GIVEN: A real data repo and plink2 conda env
    # And a snakefile with sample call rate filter 1 -> snp call rate filter 1
    conda_envs.copy_env("plink2", tmp_path)
    data_repo = (
        RealData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .copy("production_outputs/plink_start", "plink_start")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            rule all:
                input: expand("plink_filter_call_rate_1/snps.{ext}", ext=["bed", "bim", "fam"])

            rule sample_call_rate_filter_1:
                input: expand("plink_start/samples.{ext}", ext=["bed", "bim", "fam"])
                params:
                    mind=1 - cfg.config.software_params.samp_cr_1,
                    in_prefix="plink_start/samples",
                    out_prefix="plink_filter_call_rate_1/samples",
                output: expand("plink_filter_call_rate_1/samples.{ext}", ext=["bed", "bim", "fam"])
                log: "plink_filter_call_rate_1/samples.log",
                conda: cfg.conda("plink2.yml")
                threads: 2
                resources:
                    mem=1000,
                shell:
                    "plink "
                    "--bfile {params.in_prefix} "
                    "--mind {params.mind} "
                    "--make-bed "
                    "--threads {threads} "
                    "--memory {resources.mem} "
                    "--out {params.out_prefix}"

            rule snp_call_rate_filter_1:
                input: expand("plink_filter_call_rate_1/samples.{ext}", ext=["bed", "bim", "fam"])
                params:
                    geno=1 - cfg.config.software_params.snp_cr_1,
                    in_prefix="plink_filter_call_rate_1/samples",
                    out_prefix="plink_filter_call_rate_1/snps",
                output: expand("plink_filter_call_rate_1/snps.{ext}", ext=["bed", "bim", "fam"])
                log: "plink_filter_call_rate_1/snps.log",
                conda: cfg.conda("plink2.yml")
                threads: 2
                resources:
                    mem=1000,
                shell:
                    "plink "
                    "--bfile {params.in_prefix} "
                    "--geno {params.geno} "
                    "--make-bed "
                    "--threads {threads} "
                    "--memory {resources.mem} "
                    "--out {params.out_prefix}"
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path)

    # THEN: The file after sample -> snp filtering should match the file from the legacy workflow.
    obs_bed = tmp_path / "plink_filter_call_rate_1/snps.bed"
    exp_bed = data_repo / "production_outputs/plink_filter_call_rate_1/samples.bed"
    assert file_hashes_equal(obs_bed, exp_bed)


def parse_call_rate_log(file_name: Path) -> Tuple[int, int]:
    """Parse a plink log to find out how many snps/samples passed filtering."""
    m = re.findall(
        r".*\n(\d+) variants and (\d+) people pass filters and QC.*", file_name.read_text()
    )
    variants, samples = m[0]
    return int(variants), int(samples)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_call_rate_filter_1(tmp_path, conda_envs):
    """Check call rate filter 1. """
    # GIVEN: A real data repo and plink2 conda env
    conda_envs.copy_env("plink2", tmp_path)
    data_repo = (
        RealData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .copy("production_outputs/plink_start", "plink_start")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("call_rate_filters.smk")

            rule all:
                input: expand("plink_filter_call_rate_1/samples.{ext}", ext=["bed", "bim", "fam"])
            """
        )
    )

    # WHEN: I run snakemake
    snake_log = run_snakemake(tmp_path)

    # THEN:
    # I should remove the 3 temporary filtered SNP files.
    assert snake_log.count("Removing temporary output") == 3

    # and I should have the filtered sample files.
    outputs = [
        tmp_path / "plink_filter_call_rate_1/samples.bed",
        tmp_path / "plink_filter_call_rate_1/samples.bim",
        tmp_path / "plink_filter_call_rate_1/samples.fam",
    ]
    for file_name in outputs:
        assert file_name.exists()

    # the number of variants and samples filtered should be close to production outputs
    obs_var, obs_samples = parse_call_rate_log(tmp_path / "plink_filter_call_rate_1/samples.log")
    exp_var, exp_samples = parse_call_rate_log(
        data_repo / "production_outputs/plink_filter_call_rate_1/samples.log"
    )
    assert abs(obs_var - exp_var) <= 100
    assert abs(obs_samples - exp_samples) <= 5


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_call_rate_filter_2(tmp_path, conda_envs):
    """Check call rate filter 2. """
    # GIVEN: A real data repo and plink2 conda env
    conda_envs.copy_env("plink2", tmp_path)
    data_repo = (
        RealData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .copy("production_outputs/plink_filter_call_rate_1", "plink_filter_call_rate_1")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("call_rate_filters.smk")

            rule all:
                input: expand("plink_filter_call_rate_2/samples.{ext}", ext=["bed", "bim", "fam"])
            """
        )
    )

    # WHEN: I run snakemake
    snake_log = run_snakemake(tmp_path)

    # THEN:
    # I should remove the 3 temporary filtered SNP files.
    assert snake_log.count("Removing temporary output") == 3

    # and I should have the filtered sample files.
    outputs = [
        tmp_path / "plink_filter_call_rate_2/samples.bed",
        tmp_path / "plink_filter_call_rate_2/samples.bim",
        tmp_path / "plink_filter_call_rate_2/samples.fam",
    ]
    for file_name in outputs:
        assert file_name.exists()

    # the number of variants and samples filtered should be close to production outputs
    obs_var, obs_samples = parse_call_rate_log(tmp_path / "plink_filter_call_rate_2/samples.log")
    exp_var, exp_samples = parse_call_rate_log(
        data_repo / "production_outputs/plink_filter_call_rate_2/samples.log"
    )
    assert abs(obs_var - exp_var) <= 12000
    assert abs(obs_samples - exp_samples) <= 5


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
@pytest.mark.parametrize(
    "prefix", ["plink_start", "plink_filter_call_rate_1", "plink_filter_call_rate_2"]
)
def test_call_rate_stats(tmp_path, conda_envs, prefix):
    """Check call rate filter 2. """
    # GIVEN: A real data repo and plink2 conda env
    conda_envs.copy_env("plink2", tmp_path)
    data_repo = (
        RealData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .copy(f"production_outputs/{prefix}", prefix)
        .make_config()
        .make_snakefile(
            f"""
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.rules("call_rate_filters.smk")

            rule all:
                input: expand("{prefix}/samples.{{ext}}", ext=["imiss", "lmiss"])
            """
        )
    )
    # delete the stat files so snakemake will run
    # note: I have to do these gymnastics b/c the production workflow doesn't
    # use a consistent file name
    next((tmp_path / prefix).glob("*.imiss")).unlink()
    next((tmp_path / prefix).glob("*.lmiss")).unlink()

    # WHEN: I run snakemake
    run_snakemake(tmp_path)

    # THEN: The outputs should match the legacy workflow
    assert file_hashes_equal(
        tmp_path / f"{prefix}/samples.imiss",
        next((data_repo / f"production_outputs/{prefix}").glob("*.imiss")),
    )
    assert file_hashes_equal(
        tmp_path / f"{prefix}/samples.lmiss",
        next((data_repo / f"production_outputs/{prefix}").glob("*.lmiss")),
    )
