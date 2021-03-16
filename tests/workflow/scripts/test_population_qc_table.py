from pathlib import Path

import pytest

from cgr_gwas_qc.testing import chdir
from cgr_gwas_qc.testing.data import RealData
from cgr_gwas_qc.workflow.scripts.population_qc_table import (
    add_metadata,
    build_table,
    extract_files,
)


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
def test_add_metadata(qc_summary, population_level):
    with chdir(population_level):
        results = Path("population_level/results.done")
        controls = Path("population_level/controls.done")
        df = build_table(results, controls)
        df_w_metadata = add_metadata(df, qc_summary)

    assert (326, 18) == df_w_metadata.shape


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
