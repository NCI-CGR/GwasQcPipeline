import re

import pytest

from cgr_gwas_qc.testing import run_snakemake
from cgr_gwas_qc.testing.data import RealData


################################################################################
# Filter SNPs not in 1000 Genome (1KG)
################################################################################
@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_thousG_match(tmp_path, conda_envs):
    # GIVEN: Real data with ld_prune output and the plink2 conda environment
    conda_envs.copy_env("plink2", tmp_path)
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .make_config()
        .copy("production_outputs/ld_prune", "snp_filters/ld_prune")
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("snp_filters.smk")

            rule all:
                input: expand("snp_filters/plink_thousG_match/samples.{ext}", ext=["bed", "bim", "fam"])
            """
        )
    )

    # WHEN: I run snakemake to generate thousG filtered samples
    run_snakemake(tmp_path)

    # THEN: The same number of SNPs that match in the bim_filter_vcf step are
    # kept in the plink step.
    obs_summary = re.findall(
        r"\nNumber SNPs Match:\s+(\d+)\s*\n",
        (tmp_path / "snp_filters/plink_thousG_match/thousG_match.log").read_text(),
    )[0]

    obs_plink = re.findall(
        r"\n(\d+) variants and",
        (tmp_path / "snp_filters/plink_thousG_match/samples.log").read_text(),
    )[0]

    assert obs_summary == obs_plink

    # The same number of SNPs are kept in the production pipeline
    exp_plink = re.findall(
        r"\n(\d+) variants and",
        (data_store / "production_outputs/plink_thousG_match/samples.log").read_text(),
    )[0]

    assert obs_plink == exp_plink

    # The the flipped strands are identical with production pipeline
    obs_flip = tmp_path / "snp_filters/plink_thousG_match/thousG_snps_to_remove.txt"
    exp_flip = data_store / "production_outputs/ld_prune/thousG_rename.remove.txt"
    with open(obs_flip) as obs, open(exp_flip) as exp_:
        for obs_row, exp_row in zip(obs, exp_):
            assert obs_row.strip().split() == exp_row.strip().split()
