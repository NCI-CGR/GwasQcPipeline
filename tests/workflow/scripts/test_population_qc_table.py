from pathlib import Path

import pandas as pd
import pytest
from typer.testing import CliRunner

from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.testing.data import RealData
from cgr_gwas_qc.workflow.scripts.population_qc_table import (
    add_metadata,
    app,
    build_table,
    extract_files,
)

runner = CliRunner()


@pytest.mark.real_data
@pytest.fixture(scope="module")
def population_level(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("population_level")
    (
        RealData(tmp_path)
        .copy(
            "production_outputs/pca/EUR_subjects.eigenvec",
            "population_level/European/subjects_unrelated0.2_maf0.2_ld0.1_pruned.eigenvec",
        )
        .copy(
            "production_outputs/autosomal_heterozygosity/EUR_subjects_qc.het",
            "population_level/European/subjects.het",
        )
        .copy(
            "production_outputs/pca/EUR_subjects.eigenvec",
            "population_level/African/subjects_unrelated0.2_maf0.2_ld0.1_pruned.eigenvec",
        )
        .copy(
            "production_outputs/autosomal_heterozygosity/EUR_subjects_qc.het",
            "population_level/African/subjects.het",
        )
    )

    results = tmp_path / "population_level/results.done"
    results.write_text(
        "\n".join(
            [
                "population_level/European/subjects_unrelated0.2_maf0.2_ld0.1_pruned.eigenvec",
                "population_level/European/subjects_maf0.2_ld0.1_pruned.genome",
                "population_level/European/subjects.het",
                "population_level/African/subjects_unrelated0.2_maf0.2_ld0.1_pruned.eigenvec",
                "population_level/African/subjects_maf0.2_ld0.1_pruned.genome",
                "population_level/African/subjects.het",
            ]
        )
    )

    controls = tmp_path / "population_level/controls.done"
    controls.write_text(
        "population_level/European/controls_unrelated0.2_maf0.05_snps_autosome_cleaned.hwe"
        "population_level/African/controls_unrelated0.2_maf0.05_snps_autosome_cleaned.hwe"
    )

    return tmp_path


@pytest.mark.real_data
def test_build_table(population_level):
    with chdir(population_level):
        results = Path("population_level/results.done")
        controls = Path("population_level/controls.done")
        df = build_table(results, controls)

    assert (326, 16) == df.shape


@pytest.mark.real_data
def test_build_table_no_data(tmp_path):
    # GIVEN: empty results files
    results = tmp_path / "results.done"
    results.touch()
    controls = tmp_path / "controls.done"
    controls.touch()

    # WHEN: I run build table
    df = build_table(results, controls)

    # THEN: I should get an empty dataframe
    assert (0, 0) == df.shape


@pytest.mark.real_data
def test_add_metadata(sample_qc, population_level):
    with chdir(population_level):
        results = Path("population_level/results.done")
        controls = Path("population_level/controls.done")
        df = build_table(results, controls)
        df_w_metadata = add_metadata(df, sample_qc)

    assert (326, 18) == df_w_metadata.shape


@pytest.mark.real_data
def test_add_metadata_no_data(sample_qc, tmp_path):
    # GIVEN: empty results files and qc metadata
    results = tmp_path / "results.done"
    results.touch()
    controls = tmp_path / "controls.done"
    controls.touch()

    # WHEN: I build the population results table add the metadata
    df = build_table(results, controls)
    df_w_metadata = add_metadata(df, sample_qc)

    # THEN: I should get an empty dataframe
    assert (0, 0) == df_w_metadata.shape


@pytest.fixture(scope="module")
def fake_results_list(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("fake_results_list")
    filelist = tmp_path / "test.done"
    filelist.write_text(
        "\n".join(
            [
                "somepath/pop1/filename.test1",
                "somepath/pop1/filename.test2",
                "somepath/pop2/filename.test1",
                "somepath/pop2/filename.test2",
                "somepath/pop3/filename.test1",
            ]
        )
    )

    return filelist


def test_extract_files_populations(fake_results_list):
    for popfile in extract_files(fake_results_list):
        assert popfile.population in ["pop1", "pop2", "pop3"]
        assert popfile.suffix in ["test1", "test2"]


def test_extract_files_populations_no_data(tmp_path):
    # GIVEN: an empty done file
    (tmp_path / "test.done").touch()

    # WHEN: I run extract files it never enters into the loop b/c there are no
    # files.
    side_effect = 0
    for _ in extract_files(tmp_path / "test.done"):
        side_effect = 1

    # THEN: there are no loop side effects
    assert 0 == side_effect


@pytest.mark.real_data
def test_run_population_qc_table(sample_qc, population_level, tmp_path):
    # GIVEN: population level data and the qc metadata
    # WHEN: I run the script
    with chdir(population_level):
        res = runner.invoke(
            app,
            [
                sample_qc.as_posix(),
                (population_level / "population_level/results.done").as_posix(),
                (population_level / "population_level/controls.done").as_posix(),
                (tmp_path / "test.csv").as_posix(),
            ],
        )
    assert res.exit_code == 0

    # THEN: I should get the output table
    assert (326, 18) == pd.read_csv(tmp_path / "test.csv").shape


@pytest.mark.real_data
def test_run_population_qc_table_no_data(sample_qc, tmp_path):
    results = tmp_path / "results.done"
    results.touch()
    controls = tmp_path / "controls.done"
    controls.touch()

    res = runner.invoke(
        app,
        [
            sample_qc.as_posix(),
            results.as_posix(),
            controls.as_posix(),
            (tmp_path / "test.csv").as_posix(),
        ],
    )
    assert res.exit_code == 0
    assert (0, 18) == pd.read_csv(tmp_path / "test.csv").shape
