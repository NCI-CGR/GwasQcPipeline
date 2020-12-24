import re
from pathlib import Path
from typing import Tuple

import pandas as pd
import pytest

from cgr_gwas_qc.testing import make_snakefile, run_snakemake
from cgr_gwas_qc.testing.data import FakeData, RealData


@pytest.mark.workflow
def test_graf_shell_prefix(tmp_path):
    """Inject GRAF into the path.

    The GRAF binary is not currently in bioconda. This is a little hack to
    get things working. I download the graf package into a cached location
    and then inject that into the path.
    """
    make_snakefile(
        tmp_path,
        """
        from snakemake.shell import shell
        from cgr_gwas_qc.testing.data import Graf

        graf = Graf()
        shell.prefix(f"export PATH={str(graf)}:$PATH;")

        rule test:
            output:
                "snake_path",
                "graf_help"
            shell:
                'echo "${{PATH}}" > {output[0]} '
                '&& graf --help > {output[1]} '
        """,
    )

    run_snakemake(tmp_path)

    assert "cgr_gwas_qc/graf" in (tmp_path / "snake_path").read_text()
    assert "GRAF 2.4" in (tmp_path / "graf_help").read_text()


def test_graf_fingerprint_list(tmp_path):
    # GIVEN a fake config.yml
    (
        FakeData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()
            include: cfg.rules("ancestry.smk")

            rule all:
                input:
                    "ancestry/graf_fingerprints.txt"
            """
        )
    )

    # WHEN: run snakemake to get a list of graf fingerprint rsids
    run_snakemake(tmp_path)

    # THEN
    assert all(
        [
            x.startswith("rs")
            for x in (tmp_path / "ancestry/graf_fingerprints.txt").read_text().strip().split("\n")
        ]
    )


@pytest.mark.real_data
@pytest.mark.workflow
def test_extract_graf_fingerprint_markers(tmp_path, conda_envs):
    # GIVEN: Real data with level 2 call rates.
    conda_envs.copy_env("plink2", tmp_path)
    (
        RealData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .copy("production_outputs/plink_filter_call_rate_2", "plink_filter_call_rate_2")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()
            include: cfg.rules("ancestry.smk")

            rule all:
                input:
                    expand("ancestry/samples.{ext}", ext=["bed", "bim", "fam"])
            """
        )
    )

    # WHEN: run snakemake to filter out file fingerprints.
    run_snakemake(tmp_path)

    # THEN: All of the files should exist
    assert (tmp_path / "ancestry/samples.bed").exists()
    assert (tmp_path / "ancestry/samples.bim").exists()
    assert (tmp_path / "ancestry/samples.fam").exists()

    # The BIM variant IDs should all be rsids.
    assert all(
        [
            x.split()[1].startswith("rs")
            for x in (tmp_path / "ancestry/samples.bim").read_text().strip().split("\n")
        ]
    )

    # There should be 3,379 kept and reported variants in the plink log
    num_variants_keep = int(
        re.findall(
            r"\n(\d+) variants and .* pass filters", (tmp_path / "ancestry/samples.log").read_text()
        )[0]
    )
    assert num_variants_keep == 3379


@pytest.mark.real_data
@pytest.fixture(scope="session")
def graf_module(tmp_path_factory, conda_envs) -> Tuple[RealData, Path]:
    """Run the GRAF Module

    This step is really slow (~6 min).
    """
    tmp_path = tmp_path_factory.mktemp("slow")
    conda_envs.copy_env("graf_perl", tmp_path)
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .copy("production_outputs/plink_filter_call_rate_2", "plink_filter_call_rate_2")
        .make_config(workflow_params={"subject_id_to_use": "PI_Subject_ID"})
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()
            include: cfg.rules("ancestry.smk")

            rule all:
                input:
                    "ancestry/graf_relatedness.txt",
                    "ancestry/graf_relatedness.png",
                    "ancestry/graf_pop.txt",
                    "ancestry/graf_ancestry_calls.txt",
            """
        )
    )

    run_snakemake(tmp_path)
    assert (tmp_path / "ancestry/samples.fpg").exists()
    return data_store, tmp_path


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_graf_relatedness(graf_module):
    data_store, tmp_path = graf_module
    # GIVEN: Real data with GRAF module outputs.
    # THEN: Relatedness outputs should exist
    assert (tmp_path / "ancestry/graf_relatedness.txt").exists()
    assert (tmp_path / "ancestry/graf_relatedness.png").exists()

    # Samples marked as ID should also be marked as dups by plink ibd
    obs_graf = pd.read_csv(tmp_path / "ancestry/graf_relatedness.txt", sep="\t", comment="#")
    obs = {
        tuple(sorted(x))
        for x in obs_graf.query("`geno relation` == 'ID'")[["sample1", "sample2"]].itertuples(
            index=False
        )
    }

    exp_plink = pd.read_csv(data_store / "production_outputs/concordance/KnownReplicates.csv")
    exp_ = {
        tuple(sorted(x))
        for x in exp_plink.dropna()[["Sample_ID1", "Sample_ID2"]].itertuples(index=False)
    }
    assert obs == exp_


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_graf_ancestry(graf_module):
    data_store, tmp_path = graf_module
    # GIVEN: Real data with GRAF module outputs.
    # THEN: All of the files should exist
    assert (tmp_path / "ancestry/graf_pop.txt").exists()
    assert (tmp_path / "ancestry/graf_ancestry_calls.txt").exists()

    # Ancestry calls should match SNPweights
    sample2subject = (
        pd.read_csv(tmp_path / "ancestry/ssm.txt", sep="\t").set_index("Sample_ID").squeeze()
    )  # Mapping of Sample_ID to Subject_ID

    obs_ancestry_by_Subject_ID = (
        pd.read_csv(tmp_path / "ancestry/graf_ancestry_calls.txt", sep="\t")
        .set_index("Subject")  # Values of Subject are really Sample_IDs
        .rename(sample2subject)  # Convert Subject (i.e., Sample_IDs) to Subject_IDs
        .rename_axis("Subject_ID")
        .rename({"Computed population": "Ancestry"}, axis=1)
        .Ancestry.replace(  # Subject_ID to Ancestry mapping
            {
                "European": "EUR",
                "African": "AFR",
                "East Asian": "ASN",
                "African American": "ADMIXED_AFR",
                "Hispanic1": "ADMIXED_AFR",
            }
        )  # Update Ancestries to match SNPweight names
        .sort_index()
        .to_frame()
        .reset_index()
        .drop_duplicates()
    )

    exp_ancestry_by_Subject_ID = (
        pd.read_csv(data_store / "production_outputs/ancestry/subjects.snpweights.csv")
        .set_index("ID")
        .rename_axis("Subject_ID")
        .Ancestry.sort_index()
        .to_frame()
        .reset_index()
        .drop_duplicates()
    )

    df = obs_ancestry_by_Subject_ID.merge(
        exp_ancestry_by_Subject_ID, on="Subject_ID", suffixes=["_obs", "_exp"]
    )

    # TODO: Why are there these 16 differences.
    assert (df.Ancestry_obs == df.Ancestry_exp).sum() == 168
    assert (df.Ancestry_obs != df.Ancestry_exp).sum() == 16
