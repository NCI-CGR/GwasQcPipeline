import pytest

from cgr_gwas_qc.testing import run_snakemake
from cgr_gwas_qc.testing.comparison import file_hashes_equal
from cgr_gwas_qc.testing.data import RealData


################################################################################
# Filters
################################################################################
@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_call_rate_filters(tmp_path, conda_envs):
    """Check call rate filters."""
    # GIVEN: A real data repo and plink2 conda env
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy("legacy_outputs/plink_start/samples.bed", "samples.bed")
        .copy("legacy_outputs/plink_start/samples.bim", "samples.bim")
        .copy("legacy_outputs/plink_start/samples.fam", "samples.fam")
        .make_config()
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            rule all:
                input: "samples.sample_cr_filtered.snp_cr_filtered.bed"

            module plink:
                snakefile: cfg.modules("plink")

            use rule sample_call_rate_filter from plink with:
                params:
                    mind=1 - 0.8,
                    out_prefix="{prefix}.sample_cr_filtered",

            use rule snp_call_rate_filter from plink with:
                params:
                    geno=1 - 0.8,
                    out_prefix="{prefix}.snp_cr_filtered",
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path, keep_temp=True)

    # THEN: The final file after filtering should match the file from the legacy workflow.
    legacy = data_cache / "legacy_outputs/plink_filter_call_rate_1/samples.bed"
    snake = tmp_path / "samples.sample_cr_filtered.snp_cr_filtered.bed"
    assert file_hashes_equal(legacy, snake)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_maf_and_ld_prune(tmp_path, conda_envs):
    # GIVEN: A real data repository and a plink2 conda environmet
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.bed", "samples.bed")
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.bim", "samples.bim")
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.fam", "samples.fam")
        .make_config()
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            rule all:
                input: "samples_maf0.2_ld0.1_pruned.bed"

            module plink:
                snakefile: cfg.modules("plink")

            use rule maf_filter, ld_filter, ld from plink
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path, keep_temp=True)

    # THEN:
    legacy = data_cache / "legacy_outputs/ld_prune/samples.bed"
    snake = tmp_path / "samples_maf0.2_ld0.1_pruned.bed"
    assert file_hashes_equal(legacy, snake)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_maf_snp_autosome_filters(tmp_path, conda_envs):
    """Re-create the filtering done to generate the HWE subjects.

    Filter criteria include:

        - Keep Control IDs
        - Remove Related IDs
        - MAF 0.05
        - SNPs only
        - Autosome only
    """
    # GIVEN: A real data repository and a plink2 conda environmet
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy("legacy_outputs/split_by_pop/EUR_subjects.bed", "samples.bed")
        .copy("legacy_outputs/split_by_pop/EUR_subjects.bim", "samples.bim")
        .copy("legacy_outputs/split_by_pop/EUR_subjects.fam", "samples.fam")
        .copy("legacy_outputs/HWP/EUR_controls.txt", "samples_to_keep.txt")
        .copy(
            "legacy_outputs/remove_related/subjects_to_remove.txt", "samples.filtered_to_remove.txt"
        )
        .make_config()
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            ruleorder: keep_ids > remove_ids

            rule all:
                input: "samples.filtered.filtered_maf0.05_snps_autosome.keep.bed"

            module plink:
                snakefile: cfg.modules("plink")

            use rule maf_filter, snps_only_filter, autosome_only_filter, keep_ids, remove_ids, keep_bfile from plink
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path, keep_temp=True)

    # THEN:
    legacy = data_cache / "legacy_outputs/HWP/EUR_subjects.bed"
    snake = tmp_path / "samples.filtered.filtered_maf0.05_snps_autosome.keep.bed"
    assert file_hashes_equal(legacy, snake)


################################################################################
# Converters
################################################################################
@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_rename_ids(tmp_path, conda_envs):
    # GIVEN: A real data repository and a plink2 conda environmet
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy("legacy_outputs/subject_level/samples.bed", "samples.bed")
        .copy("legacy_outputs/subject_level/samples.bim", "samples.bim")
        .copy("legacy_outputs/subject_level/samples.fam", "samples.fam")
        .copy("legacy_outputs/subject_level/renameSampToSub.txt", "samples.id_map.txt")
        .make_config()
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            rule all:
                input: "samples.renamed.bed"

            module plink:
                snakefile: cfg.modules("plink")

            use rule rename_ids from plink
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path, keep_temp=True)

    # THEN:
    legacy = data_cache / "legacy_outputs/subject_level/subjects.bed"
    snake = tmp_path / "samples.renamed.bed"
    assert file_hashes_equal(legacy, snake)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_bed_to_ped(tmp_path, conda_envs):
    # GIVEN: A real data repository and a plink2 conda environmet
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy("legacy_outputs/split_by_pop/EUR_subjects.bed", "samples.bed")
        .copy("legacy_outputs/split_by_pop/EUR_subjects.bim", "samples.bim")
        .copy("legacy_outputs/split_by_pop/EUR_subjects.fam", "samples.fam")
        .copy("legacy_outputs/pca/EUR_ldPruneList.prune.in", "samples_ld0.1.prune.in")
        .copy(
            "legacy_outputs/remove_related/subjects_to_remove.txt",
            "samples_ld0.1_pruned_to_remove.txt",
        )
        .make_config()
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            rule all:
                input: "samples_ld0.1_pruned.filtered.ped"

            module plink:
                snakefile: cfg.modules("plink")

            use rule remove_ids, ld_filter, bed_to_ped from plink
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path, keep_temp=True)

    # THEN:
    legacy = data_cache / "legacy_outputs/pca/EUR_subjects_ld_prune.ped"
    snake = tmp_path / "samples_ld0.1_pruned.filtered.ped"
    assert file_hashes_equal(legacy, snake)


################################################################################
# Stats
################################################################################
@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_miss(tmp_path, conda_envs):
    # GIVEN: A real data repository and a plink2 conda environmet
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.bed", "samples.bed")
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.bim", "samples.bim")
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.fam", "samples.fam")
        .make_config()
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            rule all:
                input:
                    "samples.imiss",
                    "samples.lmiss"

            module plink:
                snakefile: cfg.modules("plink")

            use rule miss from plink
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path, keep_temp=True)

    # THEN:
    legacy = data_cache / "legacy_outputs/plink_filter_call_rate_2/samples_filter2.imiss"
    snake = tmp_path / "samples.imiss"
    assert file_hashes_equal(legacy, snake)

    legacy = data_cache / "legacy_outputs/plink_filter_call_rate_2/samples_filter2.lmiss"
    snake = tmp_path / "samples.lmiss"
    assert file_hashes_equal(legacy, snake)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_sexcheck(tmp_path, conda_envs):
    # GIVEN: A real data repository and a plink2 conda environmet
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.bed", "samples.bed")
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.bim", "samples.bim")
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.fam", "samples.fam")
        .make_config()
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            rule all:
                input:
                    "samples.sexcheck",

            module plink:
                snakefile: cfg.modules("plink")

            use rule sexcheck from plink
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path, keep_temp=True)

    # THEN:
    legacy = data_cache / "legacy_outputs/plink_filter_call_rate_2/samples_filter2.sexcheck"
    snake = tmp_path / "samples.sexcheck"
    assert file_hashes_equal(legacy, snake)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_frq(tmp_path, conda_envs):
    # GIVEN: A real data repository and a plink2 conda environmet
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.bed", "samples.bed")
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.bim", "samples.bim")
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.fam", "samples.fam")
        .make_config()
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            rule all:
                input:
                    "samples.frq",

            module plink:
                snakefile: cfg.modules("plink")

            use rule frq from plink
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path, keep_temp=True)

    # THEN:
    legacy = data_cache / "legacy_outputs/plink_filter_call_rate_2/samples_filter2.frq"
    snake = tmp_path / "samples.frq"
    assert file_hashes_equal(legacy, snake)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_hwe(tmp_path, conda_envs):
    # GIVEN: A real data repository and a plink2 conda environmet
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.bed", "samples.bed")
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.bim", "samples.bim")
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.fam", "samples.fam")
        .make_config()
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            rule all:
                input:
                    "samples.hwe",

            module plink:
                snakefile: cfg.modules("plink")

            use rule hwe from plink
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path, keep_temp=True)

    # THEN:
    legacy = data_cache / "legacy_outputs/plink_filter_call_rate_2/samples_filter2.hwe"
    snake = tmp_path / "samples.hwe"
    assert file_hashes_equal(legacy, snake)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_het(tmp_path, conda_envs):
    # GIVEN: A real data repository and a plink2 conda environmet
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.bed", "samples.bed")
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.bim", "samples.bim")
        .copy("legacy_outputs/plink_filter_call_rate_2/samples.fam", "samples.fam")
        .make_config()
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            rule all:
                input:
                    "samples.het",

            module plink:
                snakefile: cfg.modules("plink")

            use rule het from plink
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path, keep_temp=True)

    # THEN:
    legacy = data_cache / "legacy_outputs/plink_filter_call_rate_2/samples_filter2.het"
    snake = tmp_path / "samples.het"
    assert file_hashes_equal(legacy, snake)


@pytest.mark.real_data
@pytest.mark.regression
@pytest.mark.workflow
def test_genome(tmp_path, conda_envs):
    # GIVEN: A real data repository and a plink2 conda environmet
    conda_envs.copy_env("plink2", tmp_path)
    data_cache = (
        RealData(tmp_path)
        .copy("legacy_outputs/ld_prune/samples.bed", "samples.bed")
        .copy("legacy_outputs/ld_prune/samples.bim", "samples.bim")
        .copy("legacy_outputs/ld_prune/samples.fam", "samples.fam")
        .make_config()
        .make_cgr_sample_sheet()
        .make_snakefile(
            """
            from cgr_gwas_qc import load_config

            cfg = load_config()

            rule all:
                input:
                    "samples.genome",

            module plink:
                snakefile: cfg.modules("plink")

            use rule genome from plink
            """
        )
    )

    # WHEN: I run snakemake
    run_snakemake(tmp_path, keep_temp=True)

    # THEN:
    legacy = data_cache / "legacy_outputs/ibd/samples.genome"
    snake = tmp_path / "samples.genome"
    assert file_hashes_equal(legacy, snake)
