import re

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
    assert "GRAF 2.4" in (tmp_path / "snake_path").read_text()


def test_graf_fingerprint_list(tmp_path):
    # GIVEN a fake config.yml
    (
        FakeData(tmp_path)
        .add_sample_sheet()
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
@pytest.mark.workflow
def test_run_graf(tmp_path, conda_envs):
    # GIVEN: Real data with level 2 call rates.
    conda_envs.copy_env("graf_perl", tmp_path)
    (
        RealData(tmp_path)
        .add_sample_sheet()
        .add_reference_files(copy=False)
        # Note: Here I am using the CR2 files directly instead of running them
        # through rsid conversion. I am doing this because `rename_to_thousG`
        # is very slow making this test run ~5 min. Don't do this in production!
        .copy("production_outputs/plink_filter_call_rate_2", "ancestry")
        .make_config()
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

    # WHEN: run snakemake to filter out file fingerprints.
    run_snakemake(tmp_path)

    # THEN: All of the files should exist
    assert (tmp_path / "ancestry/samples.fpg").exists()
    assert (tmp_path / "ancestry/graf_relatedness.txt").exists()
    assert (tmp_path / "ancestry/graf_relatedness.png").exists()
    assert (tmp_path / "ancestry/graf_pop.txt").exists()
    assert (tmp_path / "ancestry/graf_ancestry_calls.txt").exists()
