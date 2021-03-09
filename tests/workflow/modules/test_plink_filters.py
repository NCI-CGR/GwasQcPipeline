import re
from pathlib import Path

import pytest

from cgr_gwas_qc import load_config
from cgr_gwas_qc.testing import chdir, file_hashes_equal, run_snakemake
from cgr_gwas_qc.testing.data import RealData


################################################################################
# Call Rate Filters
################################################################################
@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_old_call_rate_filter_order(tmp_path, conda_envs):
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
    data_cache = (
        RealData(tmp_path)
        .copy_sample_sheet()
        .copy("production_outputs/plink_start/samples.bed", "sample_level/samples.bed")
        .copy("production_outputs/plink_start/samples.bim", "sample_level/samples.bim")
        .copy("production_outputs/plink_start/samples.fam", "sample_level/samples.fam")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            rule all:
                input:
                    "sample_level/call_rate_2/snps.bed"

            def _sample_call_rate_filter(wildcards):
                if wildcards.cr == "2":
                    return {
                        "bed": "{prefix}/call_rate_1/snps.bed",
                        "bim": "{prefix}/call_rate_1/snps.bim",
                        "fam": "{prefix}/call_rate_1/snps.fam",
                    }

                return {
                    "bed": "{prefix}/samples.bed",
                    "bim": "{prefix}/samples.bim",
                    "fam": "{prefix}/samples.fam",
                }

            rule sample_call_rate_filter:
                input:
                    unpack(_sample_call_rate_filter)
                params:
                    mind=lambda wc: 1 - float(cfg.config.software_params.dict().get(f"samp_cr_{wc.cr}")),
                    out_prefix="{prefix}/call_rate_{cr}/samples"
                output:
                    bed="{prefix}/call_rate_{cr}/samples.bed",
                    bim="{prefix}/call_rate_{cr}/samples.bim",
                    fam="{prefix}/call_rate_{cr}/samples.fam",
                wildcard_constraints:
                    samp_cr="[1, 2]"
                log:
                    "{prefix}/call_rate_{cr}/samples.log",
                conda:
                    cfg.conda("plink2.yml")
                shell:
                    "plink "
                    "--bed {input.bed} "
                    "--bim {input.bim} "
                    "--fam {input.fam} "
                    "--mind {params.mind} "
                    "--make-bed "
                    "--out {params.out_prefix}"

            rule snp_call_rate_filter:
                input:
                    bed="{prefix}/call_rate_{cr}/samples.bed",
                    bim="{prefix}/call_rate_{cr}/samples.bim",
                    fam="{prefix}/call_rate_{cr}/samples.fam",
                params:
                    geno=lambda wc: 1 - float(cfg.config.software_params.dict().get(f"snp_cr_{wc.cr}")),
                    out_prefix="{prefix}/call_rate_{cr}/snps",
                output:
                    bed="{prefix}/call_rate_{cr}/snps.bed",
                    bim="{prefix}/call_rate_{cr}/snps.bim",
                    fam="{prefix}/call_rate_{cr}/snps.fam",
                wildcard_constraints:
                    snp_cr="[1, 2]"
                log:
                    "{prefix}/call_rate_{cr}/snps.log",
                conda:
                    cfg.conda("plink2.yml")
                shell:
                    "plink "
                    "--bed {input.bed} "
                    "--bim {input.bim} "
                    "--fam {input.fam} "
                    "--geno {params.geno} "
                    "--make-bed "
                    "--out {params.out_prefix}"
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: The file after sample -> snp filtering should match the file from the legacy workflow.
    obs_cr1 = tmp_path / "sample_level/call_rate_1/snps.bed"
    exp_cr1 = data_cache / "production_outputs/plink_filter_call_rate_1/samples.bed"
    assert file_hashes_equal(obs_cr1, exp_cr1)

    obs_cr2 = tmp_path / "sample_level/call_rate_2/snps.bed"
    exp_cr2 = data_cache / "production_outputs/plink_filter_call_rate_2/samples.bed"
    assert file_hashes_equal(obs_cr2, exp_cr2)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_call_rate_filters(tmp_path, conda_envs):
    """Check call rate filters."""
    # GIVEN: A real data repo and plink2 conda env
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy_sample_sheet()
        .copy("production_outputs/plink_start/samples.bed", "sample_level/samples.bed")
        .copy("production_outputs/plink_start/samples.bim", "sample_level/samples.bim")
        .copy("production_outputs/plink_start/samples.fam", "sample_level/samples.fam")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("plink_filters.smk")

            rule all:
                input: expand("sample_level/call_rate_2/samples.{ext}", ext=["bed", "bim", "fam"])
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path)

    # THEN:
    # the number of variants and samples filtered should be close to production outputs
    def parse_call_rate_log(file_name: Path):
        """Parse a plink log to find out how many snps/samples passed filtering."""
        m = re.findall(
            r".*\n(\d+) variants and (\d+) people pass filters and QC.*", file_name.read_text()
        )
        variants, samples = m[0]
        return int(variants), int(samples)

    obs_var, obs_samples = parse_call_rate_log(tmp_path / "sample_level/call_rate_2/samples.log")
    exp_var, exp_samples = parse_call_rate_log(
        data_cache / "production_outputs/plink_filter_call_rate_2/samples.log"
    )
    assert abs(obs_var - exp_var) <= 12000
    assert abs(obs_samples - exp_samples) <= 5


################################################################################
# MAF and LD Filters
################################################################################
@pytest.fixture(scope="session")
def maf_and_ld_outputs(tmp_path_factory, conda_envs):
    tmp_path = tmp_path_factory.mktemp("maf_and_ld_outputs")

    # GIVEN: A real data repository and a plink2 conda environmet
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy_sample_sheet()
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples.bed",
            "sample_level/call_rate_2/samples.bed",
        )
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples.bim",
            "sample_level/call_rate_2/samples.bim",
        )
        .copy(
            "production_outputs/plink_filter_call_rate_2/samples.fam",
            "sample_level/call_rate_2/samples.fam",
        )
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("plink_filters.smk")

            rule all:
                input:
                    "sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.bed".format(
                        maf=cfg.config.software_params.maf_for_ibd,
                        ld=cfg.config.software_params.ld_prune_r2
                    )
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path, keep_temp=True)

    return data_cache, tmp_path


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_maf_filter(maf_and_ld_outputs):
    # GIVEN: Outputs from the maf_filter rule
    data_cache, tmp_path = maf_and_ld_outputs

    # THEN: The number of MAF filtered snps should be the same as the legacy run
    # Note: the legacy workflow does multiple filtering steps at the same time,
    # so I need to parse the log and just pull out the MAF filtering.
    def parse_plink_log_for_maf_counts(file_name) -> int:
        m = re.findall(
            r"\n(\d+) variants removed due to minor allele threshold", file_name.read_text()
        )
        return int(m[0])

    with chdir(tmp_path):
        cfg = load_config()
        maf = cfg.config.software_params.maf_for_ibd

    obs_maf_count = parse_plink_log_for_maf_counts(
        tmp_path / f"sample_level/call_rate_2/samples_maf{maf}.log"
    )
    exp_maf_count = parse_plink_log_for_maf_counts(
        data_cache / "production_outputs/ld_prune/ldPruneList.log"
    )
    assert obs_maf_count == exp_maf_count


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_approx_ld(maf_and_ld_outputs):
    # GIVEN: A real data repository and a plink2 conda environmet
    data_cache, tmp_path = maf_and_ld_outputs

    # THEN: I expect the number of pruned snps to be the same.
    def parse_plink_log_for_ld_pruning(file_name):
        m = re.findall(
            r"Pruning complete.\s+(\d+) of (\d+) variants removed.", file_name.read_text()
        )
        pruned, kept = m[0]
        return int(pruned), int(kept)

    with chdir(tmp_path):
        cfg = load_config()
        maf, ld = cfg.config.software_params.maf_for_ibd, cfg.config.software_params.ld_prune_r2

    obs_pruned, _ = parse_plink_log_for_ld_pruning(
        tmp_path / f"sample_level/call_rate_2/samples_maf{maf}_ld{ld}.log"
    )
    exp_pruned, _ = parse_plink_log_for_ld_pruning(
        data_cache / "production_outputs/ld_prune/ldPruneList.log"
    )
    assert obs_pruned == exp_pruned


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_ld_prune(maf_and_ld_outputs):
    # GIVEN:
    data_cache, tmp_path = maf_and_ld_outputs

    # THEN:
    with chdir(tmp_path):
        cfg = load_config()
        maf, ld = cfg.config.software_params.maf_for_ibd, cfg.config.software_params.ld_prune_r2

    obs_bed = tmp_path / f"sample_level/call_rate_2/samples_maf{maf}_ld{ld}_pruned.bed"
    exp_bed = data_cache / "production_outputs/ld_prune/samples.bed"
    assert file_hashes_equal(obs_bed, exp_bed)
