import re
from typing import Tuple

import pandas as pd
import pytest

from cgr_gwas_qc.testing import run_snakemake
from cgr_gwas_qc.testing.data import FakeData, RealData


@pytest.mark.workflow
def test_Sample_IDs_above_contam_threshold(tmp_path):
    # GIVEN: Fake data with the GTC entry point
    (
        FakeData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("common.smk")
            include: cfg.modules("sample_contamination.smk")

            rule all:
                input:
                    "sample_contamination/contaminated_samples.txt"
            """
        )
    )

    # Add a fake agg_contamination_test.csv table with one sample with a %Mix
    # over the contamination threshold (default 0.2).
    (tmp_path / "sample_contamination").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            ["T0001", 0.05, 0.0, 0.0],
            ["T0002", 0.1, 0.0, 0.0],
            ["T0003", 0.2, 0.0, 0.0],
            ["T0004", 0.25, 0.0, 0.0],
        ],
        columns=["Sample_ID", "%Mix", "LLK", "LLK0"],
    ).to_csv(tmp_path / "sample_contamination/agg_contamination_test.csv", index=False)

    # WHEN: I run snakemake to create a list of contaminated samples
    run_snakemake(tmp_path)

    # THEN: I expect a single sample to be flagged as contaminated.
    obs_file = tmp_path / "sample_contamination/contaminated_samples.txt"
    obs = obs_file.read_text().strip().split("\n")
    assert len(obs) == 1


def parse_plink_contamination_filtering(file_name) -> Tuple[int, int]:
    m = re.findall(
        r"\n(\d+) people.*\n--remove: (\d+) people remaining", file_name.read_text(), re.DOTALL
    )
    n_start, n_remaining = m[0]
    return int(n_start), int(n_remaining)


@pytest.mark.real_data
@pytest.mark.workflow
def test_remove_contaminated_samples(tmp_path, conda_envs):
    # GIVEN: Real data with the GTC entry point, plink files after the second
    # call rate filter, and the plink conda environment
    conda_envs.copy_env("plink2", tmp_path)
    data_store = (
        RealData(tmp_path)
        .add_sample_sheet(full_sample_sheet=False)
        .add_reference_files(copy=False)
        .add_user_files(entry_point="gtc", copy=False)
        .copy("production_outputs/plink_filter_call_rate_2", "plink_filter_call_rate_2")
        .make_config()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            include: cfg.modules("common.smk")
            include: cfg.modules("sample_contamination.smk")

            rule all:
                input:
                    expand("sample_contamination/not_contaminated/samples.{ext}", ext=["bed", "bim",  "fam"])
            """
        )
    )

    # and a contaminated sample list with 1 sample id.
    contam_sample = data_store.ss.data.Sample_ID[0]
    (tmp_path / "sample_contamination").mkdir(parents=True, exist_ok=True)
    (tmp_path / "sample_contamination/contaminated_samples.txt").write_text(
        "{} {}\n".format(contam_sample, contam_sample)
    )

    # WHEN: I run snakemake to remove contaminated samples
    run_snakemake(tmp_path)

    # THEN: I expect the files will be created
    outputs = [
        tmp_path / "sample_contamination/not_contaminated/samples.bed",
        tmp_path / "sample_contamination/not_contaminated/samples.bim",
        tmp_path / "sample_contamination/not_contaminated/samples.fam",
    ]
    assert all(file_name.exists() for file_name in outputs)

    # and the plink log shows 1 sample was filtered
    n_start, n_remaining = parse_plink_contamination_filtering(
        tmp_path / "sample_contamination/not_contaminated/samples.log"
    )
    assert (n_start - n_remaining) == 1
