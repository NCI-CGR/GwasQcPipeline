#!/usr/bin/env python
"""Helper script to compare results with the legacy workflow. """
from collections import namedtuple
from functools import partial
from pathlib import Path
from textwrap import dedent
from typing import Any, Set

import pandas as pd
import typer
from pandas.testing import assert_frame_equal, assert_series_equal
from scipy.stats import pearsonr
from typer.colors import BLUE, BRIGHT_CYAN, GREEN, MAGENTA, RED, WHITE, YELLOW

from cgr_gwas_qc import load_config, yaml
from cgr_gwas_qc.models.config import Config
from cgr_gwas_qc.testing import chdir, comparison
from cgr_gwas_qc.workflow.scripts import sample_qc_table, subject_qc_table

app = typer.Typer(add_completion=False)

Comparison = namedtuple("Comparison", "legacy,current,message,func")


@app.command()
def main(
    legacy_dir: Path = typer.Argument(..., help="Path to a legacy run using these same data."),
    qc_dir: Path = typer.Argument(
        ".",
        help="Path to a run of the CGR GwasQcPipeline. Defaults to the current working directory.",
    ),
    ignored_config: bool = typer.Option(
        False, help="True to show which legacy config options are ignored for comparisons."
    ),
    config_only: bool = typer.Option(False, help="Only compare the config."),
    mix_atol: float = typer.Option(
        0.001, help="The absolute tolerance level when comparing verifyIDintensities %Mix column."
    ),
    llk_atol: float = typer.Option(
        500.0, help="The absolute tolerance level when comparing verifyIDintensities LLK column."
    ),
    llk0_atol: float = typer.Option(
        500.0, help="The absolute tolerance level when comparing verifyIDintensities LLK0 column."
    ),
):
    """Compare results with a legacy run.

    Runs a comparison with the legacy QC workflow.
    """
    with chdir(qc_dir):
        cfg = load_config()
        typer.secho("# Check Configuration", fg=BRIGHT_CYAN)
        compare_config(cfg.config, legacy_dir, ignored_config)

        if config_only:
            typer.secho(
                "Only ran config comparison. Remove --config-only to run full comparison.",
                fg=YELLOW,
            )
            raise typer.Exit()

        typer.secho("\n# Check Sample Level Results", fg=BRIGHT_CYAN)
        compare_initial_plink(legacy_dir)
        compare_cr1_plink(legacy_dir)
        compare_cr2_plink(legacy_dir)
        compare_concordance(
            legacy_dir,
            cfg.config.software_params.maf_for_ibd,
            cfg.config.software_params.ld_prune_r2,
        )
        compare_contamination(legacy_dir, mix_atol, llk_atol, llk0_atol)
        compare_ancestry(legacy_dir)
        compare_sample_analytic_exclusions(legacy_dir)

        typer.secho("\n# Check Subject Level Results", fg=BRIGHT_CYAN)
        compare_selected_subjects(legacy_dir)
        compare_subject_analytic_exclusions(legacy_dir)

        typer.secho("\n# Compared Table 4a", fg=BRIGHT_CYAN)
        compare_4a(legacy_dir)

        typer.secho("\n# Compared Table 4b", fg=BRIGHT_CYAN)
        compare_4b(legacy_dir)


def _compare_config_values(legacy, current, name, problems):
    try:
        assert legacy == current
    except AssertionError:
        problems.add((legacy, current, name))
        typer.secho(f"The values for {name} did not match: {legacy} vs {current}", fg=RED)


def _ignored_legacy_options(legacy_config):
    """Discontinued Options.

    These options don't really have a similar value in the current workflow. I
    just want to make it clear that I am ignoring them.
    """
    _ignored = [
        (k := "plink_genotype_file", legacy_config[k]),
        (k := "lims_output_dir", legacy_config[k]),
        (k := "ibd_pi_hat_cutoff", legacy_config[k]),
        (k := "adpc_file", legacy_config[k]),
        (k := "gtc_dir", legacy_config[k]),
        (k := "remove_sex_discordant", legacy_config[k]),
        (k := "remove_unexpected_rep", legacy_config[k]),
        (k := "doc_template", legacy_config[k]),
        (k := "start_time", legacy_config[k]),
    ]
    for opt, value in _ignored:
        typer.secho(f"Ignored legacy option: {opt} [{value}]", fg=YELLOW)


def compare_config(config: Config, legacy_dir: Path, ignored_config: bool):
    legacy_config = yaml.load(legacy_dir / "config.yaml")
    options = [
        (
            Path(legacy_config["sample_sheet"]).resolve(),
            Path(config.sample_sheet).resolve(),
            "Sample Sheet",
        ),
        (
            Path(legacy_config["illumina_manifest_file"]).resolve(),
            config.reference_files.illumina_manifest_file.resolve(),
            "Illumina BPM File",
        ),
        (
            legacy_config["subject_id_to_use"],
            config.workflow_params.subject_id_column,
            "Subject ID Column",
        ),
        (
            legacy_config["expected_sex_col_name"],
            config.workflow_params.expected_sex_column,
            "Expected Sex Column",
        ),
        (legacy_config.get("num_samples", "NA"), config.num_samples, "Number of Samples"),
        (legacy_config.get("numSNPs", "NA"), config.num_snps, "Number of SNPs"),
        (legacy_config["snp_cr_1"], config.software_params.snp_call_rate_1, "SNP Call Rate 1"),
        (
            legacy_config["samp_cr_1"],
            config.software_params.sample_call_rate_1,
            "Sample Call Rate 1",
        ),
        (legacy_config["snp_cr_2"], config.software_params.snp_call_rate_2, "SNP Call Rate 2"),
        (
            legacy_config["samp_cr_2"],
            config.software_params.sample_call_rate_2,
            "Sample Call Rate 2",
        ),
        (legacy_config["maf_for_ibd"], config.software_params.maf_for_ibd, "MAF for IBD"),
        (legacy_config["ld_prune_r2"], config.software_params.ld_prune_r2, "LD Prune R^2"),
        (
            legacy_config["contam_threshold"],
            config.software_params.contam_threshold,
            "Contamination Threshold",
        ),
        (
            legacy_config["contam_pop"],
            config.software_params.contam_population,
            "B Allele Frequency Population (contamination check)",
        ),
        (
            legacy_config["strand"].lower(),
            config.software_params.strand,
            "Strand",
        ),
        (
            legacy_config["pi_hat_threshold"],
            config.software_params.pi_hat_threshold,
            "PI_HAT Threshold",
        ),
        (
            legacy_config["dup_concordance_cutoff"],
            config.software_params.dup_concordance_cutoff,
            "Concordance Threshold",
        ),
        (
            legacy_config["autosomal_het_thresh"],
            config.software_params.autosomal_het_threshold,
            "Autosomal Heterozygosity Threshold",
        ),
        (
            legacy_config["minimum_pop_subjects"],
            config.workflow_params.minimum_pop_subjects,
            "Minimum Number of Subjects Per Population",
        ),
        (
            legacy_config["control_hwp_thresh"],
            config.workflow_params.control_hwp_threshold,
            "Minimum Number of Controls Per Population",
        ),
        (
            bool(legacy_config["remove_contam"]),
            config.workflow_params.remove_contam,
            "Remove Contaminated Samples",
        ),
        (
            bool(legacy_config["remove_rep_discordant"]),
            config.workflow_params.remove_rep_discordant,
            "Remove Replicate Discordant Samples",
        ),
    ]
    problems: Set[Any] = set()
    for option in options:
        _compare_config_values(*option, problems)

    if not problems:
        typer.secho("Config Options Match", fg=GREEN)

    if ignored_config:
        _ignored_legacy_options(legacy_config)


def _file_comparison(cmp):
    try:
        match = "File Hashes Match" if cmp.func.__name__ == "file_hashes_equal" else "File Match"
        cmp.func(cmp.legacy, cmp.current)
        typer.secho(f"{cmp.message} {match}", fg=GREEN)
    except AssertionError:
        typer.secho(f"{cmp.message} did not match ({cmp.legacy} vs {cmp.current}).", fg=RED)
    except NotImplementedError:
        typer.secho(
            f"Cannot currently compare {cmp.message}",
            fg=YELLOW,
        )


def compare_initial_plink(legacy_dir: Path):
    typer.secho("\n## Comparing Initial PLINK Files", fg=MAGENTA)

    files = [
        Comparison(
            legacy_dir / "plink_start/samples.bed",
            "sample_level/samples.bed",
            "Initial BED",
            comparison.assert_plink_bed_equal,
        ),
        Comparison(
            legacy_dir / "plink_start/samples.bim",
            "sample_level/samples.bim",
            "Initial BIM",
            comparison.assert_plink_bim_equal,
        ),
        Comparison(
            legacy_dir / "plink_start/samples.fam",
            "sample_level/samples.fam",
            "Initial FAM",
            comparison.assert_plink_fam_equal,
        ),
        Comparison(
            legacy_dir / "plink_start/samples_start.imiss",
            "sample_level/samples.imiss",
            "Initial IMISS",
            comparison.file_hashes_equal,
        ),
        Comparison(
            legacy_dir / "plink_start/samples_start.lmiss",
            "sample_level/samples.lmiss",
            "Initial LMISS",
            comparison.file_hashes_equal,
        ),
    ]
    [_file_comparison(cmp) for cmp in files]


def compare_cr1_plink(legacy_dir: Path):
    typer.secho("\n## Comparing Call Rate 1 PLINK Files", fg=MAGENTA)

    files = [
        Comparison(
            legacy_dir / "plink_filter_call_rate_1/samples.bed",
            "sample_level/call_rate_1/samples.bed",
            "Call Rate 1 BED",
            comparison.assert_plink_bed_equal,
        ),
        Comparison(
            legacy_dir / "plink_filter_call_rate_1/samples.bim",
            "sample_level/call_rate_1/samples.bim",
            "Call Rate 1 BIM",
            comparison.assert_plink_bim_equal,
        ),
        Comparison(
            legacy_dir / "plink_filter_call_rate_1/samples.fam",
            "sample_level/call_rate_1/samples.fam",
            "Call Rate 1 FAM",
            comparison.assert_plink_fam_equal,
        ),
        Comparison(
            legacy_dir / "plink_filter_call_rate_1/samples_filter1.imiss",
            "sample_level/call_rate_1/samples.imiss",
            "Call Rate 1 IMISS",
            comparison.file_hashes_equal,
        ),
        Comparison(
            legacy_dir / "plink_filter_call_rate_1/samples_filter1.lmiss",
            "sample_level/call_rate_1/samples.lmiss",
            "Call Rate 1 LMISS",
            comparison.file_hashes_equal,
        ),
        Comparison(
            legacy_dir / "plink_filter_call_rate_1/samples_filter1.sexcheck",
            "sample_level/call_rate_1/samples.sexcheck",
            "Call Rate 1 SexCheck",
            comparison.file_hashes_equal,
        ),
    ]
    [_file_comparison(cmp) for cmp in files]


def compare_cr2_plink(legacy_dir: Path):
    typer.secho("\n## Comparing Call Rate 2 PLINK Files", fg=MAGENTA)

    files = [
        Comparison(
            legacy_dir / "plink_filter_call_rate_2/samples.bed",
            "sample_level/call_rate_2/samples.bed",
            "Call Rate 2 BED",
            comparison.assert_plink_bed_equal,
        ),
        Comparison(
            legacy_dir / "plink_filter_call_rate_2/samples.bim",
            "sample_level/call_rate_2/samples.bim",
            "Call Rate 2 BIM",
            comparison.assert_plink_bim_equal,
        ),
        Comparison(
            legacy_dir / "plink_filter_call_rate_2/samples.fam",
            "sample_level/call_rate_2/samples.fam",
            "Call Rate 2 FAM",
            comparison.assert_plink_fam_equal,
        ),
        Comparison(
            legacy_dir / "plink_filter_call_rate_2/samples_filter2.imiss",
            "sample_level/call_rate_2/samples.imiss",
            "Call Rate 2 IMISS",
            comparison.file_hashes_equal,
        ),
        Comparison(
            legacy_dir / "plink_filter_call_rate_2/samples_filter2.lmiss",
            "sample_level/call_rate_2/samples.lmiss",
            "Call Rate 2 LMISS",
            comparison.file_hashes_equal,
        ),
    ]
    [_file_comparison(cmp) for cmp in files]


def compare_concordance(legacy_dir: Path, maf: float, ld: float):
    typer.secho("\n## Comparing Sample Concordance", fg=MAGENTA)
    _file_comparison(
        Comparison(
            legacy_dir / "ibd/samples.genome",
            f"sample_level/call_rate_2/samples_maf{maf}_ld{ld}.genome",
            "Sample Level IBD",
            comparison.file_hashes_equal,
        )
    )


def _contamination_differences(cmp: Comparison, mix_atol, llk_atol, llk0_atol) -> str:
    legacy = pd.read_csv(cmp.legacy, index_col="ID").rename_axis("Sample_ID")
    dev = pd.read_csv(cmp.current, index_col="Sample_ID").drop(
        "is_contaminated", axis=1, errors="ignore"
    )

    # If median IDAT intensity is <6000 and call rate is low then the workflows
    # should set %Mix to NA. To indicate no difference between workflows, I am
    # setting NA to 0 when both workflows have %Mix == NA.
    both_na = legacy["%Mix"].isna() & dev["%Mix"].isna()
    legacy.loc[both_na, "%Mix"] = 0
    dev.loc[both_na, "%Mix"] = 0

    df = pd.concat(
        [
            (legacy["%Mix"] - dev["%Mix"]).abs().rename("Mix_diff"),
            (legacy.LLK - dev.LLK).abs().rename("LLK_diff"),
            (legacy.LLK0 - dev.LLK0).abs().rename("LLK0_diff"),
        ],
        axis=1,
    )

    num_missing = df["Mix_diff"].isna().sum()
    suffix = ""
    if num_missing > 0:
        suffix = dedent(
            f"""
            There are {num_missing:,} rows where %Mix was set to NA in one
            workflow but not the other. Both workflows should set %Mix to NA if
            IDAT intensity is less than the given threshold and the call rate
            was low.
            """
        )

    return (
        "\n"
        + df.query(
            "Mix_diff.isna() | Mix_diff > @mix_atol | "
            "LLK_diff.isna() | LLK_diff > @llk_atol | "
            "LLK0_diff.isna() | LLK0_diff > @llk0_atol"
        )
        .sort_values(["Mix_diff", "LLK_diff", "LLK0_diff"])
        .to_string(max_rows=20, show_dimensions=True)
        + "\n\n"
        + df.agg(["min", "max"]).to_string()
        + "\n\n"
        + suffix
    )


def compare_contamination(legacy_dir: Path, mix_atol: float, llk_atol: float, llk0_atol: float):
    typer.secho("\n## Comparing Contamination (Fuzzy Match)", fg=MAGENTA)
    typer.secho(
        (
            "Fuzzy match thresholds:\n"
            f"  - %Mix +/- {mix_atol}\n"
            f"  - LLK +/- {llk_atol}\n"
            f"  - LLK0 +/- {llk0_atol}\n"
        ),
        fg=YELLOW,
    )

    cmp = Comparison(
        legacy_dir / "all_contam/contam.csv",
        "sample_level/contamination/summary.csv",
        "verifyIDintensity Results",
        partial(
            comparison.assert_verifyIDintensity_equal,
            mix_atol=mix_atol,
            llk_atol=llk_atol,
            llk0_atol=llk0_atol,
        ),
    )

    try:
        cmp.func(cmp.legacy, cmp.current)
        typer.secho(f"{cmp.message} Fuzzy Match", fg=GREEN)
    except AssertionError:
        typer.secho(
            f"{cmp.message} did not match at current tolerance rates ({cmp.legacy} vs {cmp.current}). "
            "You may want to adjust the tolerance rates using `--mix-atol`, `--llk-atol`, or "
            "`--llk0-atol`. Below is summary of the absolute differences.",
            fg=RED,
        )
        typer.secho(_contamination_differences(cmp, mix_atol, llk_atol, llk0_atol))


def _parse_snpweights(filename):
    return (
        pd.read_csv(filename, dtype={"ID": str})
        .rename(
            {
                "ID": "Sample_ID",
                "AFR": "Pct_AFR",
                "ASN": "Pct_ASN",
                "EUR": "Pct_EUR",
            },
            axis=1,
        )
        .replace(
            {
                "EUR": "European",
                "AFR": "African",
                "ASN": "East Asian",
            }
        )
        .set_index("Sample_ID")
        .sort_index()
    )


def _parse_graf(filename):
    return (
        pd.read_csv(filename, sep="\t", dtype={"Subject": str})
        .rename(
            {
                "Subject": "Sample_ID",
                "Computed population": "Ancestry",
                "P_f (%)": "Pct_AFR",
                "P_a (%)": "Pct_ASN",
                "P_e (%)": "Pct_EUR",
            },
            axis=1,
        )
        .set_index("Sample_ID")
        .sort_index()
    )


def _cross_tabulate_ancestry(legacy, dev):
    typer.secho("\n### Ancestry Assignment Crosstabulation", fg=BLUE)
    table = pd.crosstab(legacy, dev).rename_axis("Legacy").rename_axis("Production", axis=1)
    typer.secho("Please carefully check the table for differences", fg=YELLOW)
    typer.secho(table.to_string())


def _count_differences_major_groups(legacy, dev):
    typer.secho("\n### Difference in the Number of Samples Assigned to Major Ancestry", fg=BLUE)
    major_ancestry_groups = ["African", "East Asian", "European"]
    legacy_major = legacy[legacy.isin(major_ancestry_groups)]
    dev_major = dev[dev.isin(major_ancestry_groups)]

    try:
        assert_series_equal(legacy_major, dev_major, check_names=False)
        typer.secho("Major Ancestry Match", fg=GREEN)
    except AssertionError:
        legacy_counts = legacy_major.value_counts().sort_index().rename("Legacy")
        dev_counts = dev_major.value_counts().sort_index().rename("Production")
        table = (
            pd.concat([legacy_counts, dev_counts], axis=1)
            .fillna(0)
            .astype(int)
            .assign(Difference=lambda x: (x.iloc[:, 0] - x.iloc[:, 1]).abs())
        )

        try:
            assert_series_equal(legacy_counts, dev_counts, check_names=False, rtol=0.05)
            typer.secho("Please carefully check the table for small differences", fg=YELLOW)
            color = WHITE
        except AssertionError:
            typer.secho("Please check the table for large differences", fg=RED)
            color = RED

        typer.secho(table.to_string(), fg=color)


def _corr_of_pct_ancestry(legacy, dev, name):
    R2 = pearsonr(legacy, dev)[0] ** 2
    color = GREEN if R2 >= 0.95 else RED
    typer.secho(f"{name} R^2: {R2:0.4}", fg=color)


def compare_ancestry(legacy_dir: Path):
    typer.secho("\n## Comparing Ancestry (Fuzzy Match)", fg=MAGENTA)
    legacy = _parse_snpweights(legacy_dir / "snpweights/samples.snpweights.csv")
    dev = _parse_graf("sample_level/ancestry/graf_ancestry.txt")
    _cross_tabulate_ancestry(legacy.Ancestry, dev.Ancestry)
    _count_differences_major_groups(legacy.Ancestry, dev.Ancestry)

    typer.secho("\n### Correlation of Percent Ancestry Estimates", fg=BLUE)
    _corr_of_pct_ancestry(legacy.Pct_AFR, dev.Pct_AFR, "% AFR")
    _corr_of_pct_ancestry(legacy.Pct_ASN, dev.Pct_ASN, "% ASN")
    _corr_of_pct_ancestry(legacy.Pct_EUR, dev.Pct_EUR, "% EUR")


def _legacy_sample_level_exclusions(filename):
    """Wrangle legacy QC table to get exclusions"""
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _get_reason

    exclusion_criteria = {
        "is_missing_idats": "Missing IDATs",
        "Call_Rate_1_filter": "Call Rate 1 Filtered",
        "Call_Rate_2_filter": "Call Rate 2 Filtered",
        "Contaminated": "Contamination",
        "Expected Replicate Discordance": "Replicate Discordance",
    }

    return (
        pd.read_csv(filename, dtype={"Sample_ID": str})
        .set_index("Sample_ID")
        .replace({"N": False, "Y": True})
        .assign(is_missing_idats=lambda x: ~x.IdatsInProjectDir.astype(bool))
        .astype({**{k: "boolean" for k in exclusion_criteria.keys()}})
        .assign(Call_Rate_2_filter=lambda x: x.Call_Rate_2_filter ^ x.Call_Rate_1_filter)
        .assign(Legacy_Exclusions=lambda x: _get_reason(x, exclusion_criteria))
        .Legacy_Exclusions.sort_index()
        .fillna("")
    )


def _sample_level_exclusions(filename):
    return (
        sample_qc_table.read(filename)
        .set_index("Sample_ID")
        .analytic_exclusion_reason.rename("Production_Exclusions")
        .sort_index()
        .fillna("")
    )


def compare_sample_analytic_exclusions(legacy_dir: Path):
    typer.secho("\n## Exclusion Reason By Sample", fg=MAGENTA)
    legacy = _legacy_sample_level_exclusions(legacy_dir / "all_sample_qc.csv")
    dev = _sample_level_exclusions("sample_level/sample_qc.csv")
    try:
        assert_series_equal(legacy, dev, check_names=False, check_dtype=False)
        typer.secho("Per Sample Analytic Exclusions Match", fg=GREEN)
    except AssertionError:
        typer.secho("Please carefully check this table\n", fg=YELLOW)
        mask = legacy != dev
        df = pd.concat([legacy[mask], dev[mask]], axis=1)
        typer.secho("Differences in Exclusion Reason", fg=RED)
        typer.secho(df.to_string(max_rows=20, show_dimensions=True))
        typer.secho("Summary of differences", fg=RED)
        typer.secho(df.value_counts().to_string())
        typer.secho(f"Total number of differences {mask.sum():,}", fg=RED)


def _legacy_selected_subjects(filename):
    return (
        pd.read_csv(filename, dtype={"Subject_ID": str, "Sample_ID": str})
        .set_index("Subject_ID")
        .sort_index()
        .squeeze()
        .rename("Legacy")
    )


def _selected_subjects(filename):
    df = (
        sample_qc_table.read(filename)
        .set_index("Group_By_Subject_ID")
        .rename_axis("Subject_ID")
        .query("not is_internal_control")
    )
    subjects = df.index.unique().sort_values()
    return df.query("is_subject_representative").Sample_ID.reindex(subjects).rename("Production")


def compare_selected_subjects(legacy_dir: Path):
    typer.secho("\n## Selected Subject Representative", fg=MAGENTA)
    legacy = _legacy_selected_subjects(legacy_dir / "subject_level/SampleUsedforSubject.csv")
    dev = _selected_subjects("sample_level/sample_qc.csv")
    try:
        assert_series_equal(legacy, dev, check_names=False, check_dtype=False)
        typer.secho("Selected Subjects Match", fg=GREEN)
    except AssertionError:
        typer.secho("Selected Subjects Different", fg=RED)
        table = pd.concat([legacy.fillna("Dropped"), dev.fillna("Dropped")], axis=1).query(
            "Legacy != Production"
        )
        typer.secho(table.to_string(), fg=RED)


def _legacy_subject_level_exclusions(filename):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _get_reason

    exclusion_criteria = {
        "Sex Discordant": "Sex Discordance",
        "Unexpected Replicate": "Unexpected Replicate",
    }

    return (
        pd.read_csv(filename, dtype={"Sample_ID": str})
        .set_index("Sample_ID")
        .replace({"N": False, "Y": True})
        .astype({**{k: "boolean" for k in exclusion_criteria.keys()}})
        .assign(Legacy_Exclusions=lambda x: _get_reason(x, exclusion_criteria))
        .Legacy_Exclusions.sort_index()
        .fillna("")
    )


def _subject_level_exclusions(filename):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _get_reason

    exclusion_criteria = {
        "is_sex_discordant": "Sex Discordance",
        "is_unexpected_replicate": "Unexpected Replicate",
    }

    return (
        subject_qc_table.read(filename)
        .set_index("Sample_ID")
        .assign(Production_Exclusions=lambda x: _get_reason(x, exclusion_criteria))
        .Production_Exclusions.sort_index()
        .fillna("")
    )


def compare_subject_analytic_exclusions(legacy_dir: Path):
    typer.secho("\n## Exclusion Reason By Subject", fg=MAGENTA)
    legacy = _legacy_subject_level_exclusions(legacy_dir / "all_sample_qc.csv")
    dev = _subject_level_exclusions("subject_level/subject_qc.csv")
    legacy = legacy.reindex(dev.index)  # drop samples that are  not subjects
    try:
        assert_series_equal(legacy, dev, check_names=False)
        typer.secho("Per Subject Analytic Exclusions Match", fg=GREEN)
    except AssertionError:
        typer.secho("Please carefully check this table\n", fg=YELLOW)
        mask = legacy != dev
        df = pd.concat([legacy[mask], dev[mask]], axis=1)
        typer.secho("Differences in Exclusion Reason", fg=RED)
        typer.secho(df.to_string(max_rows=20, show_dimensions=True))
        typer.secho("Summary of differences", fg=RED)
        typer.secho(df.value_counts().to_string())
        typer.secho(f"Total number of differences {mask.sum():,}", fg=RED)


def _legacy_table_4a(exclusions, remaining):
    counts = (
        pd.read_csv(exclusions)
        .set_index("Reason")
        .rename_axis("Filter Reason/Description")
        .drop(["SexDiscordant", "UnexpectedReplicates", "AutosomalHet"])
        .rename({"Total": "All Samples"}, axis=1)
    )

    totals = (
        pd.read_csv(remaining)
        .set_index("ExlusionPriorToCounts")
        .rename_axis("Filter Reason/Description")
        .reindex(["InternalQC", "SubjectLevel"])
        .rename(
            index={
                "InternalQC": "Samples Remaining for Analysis",
                "SubjectLevel": "Subjects Remaining",
            },
            columns={"Total": "All Samples"},
        )
    )

    return counts.append(totals)


def _table_4a(sample_qc_csv):
    from cgr_gwas_qc.reporting.qc_exclusions import sample_exclusion_counts

    return sample_exclusion_counts(sample_qc_csv).drop("User Excluded")


def compare_4a(legacy_dir: Path):
    legacy = _legacy_table_4a(
        legacy_dir / "counts/exclusion_counts.csv", legacy_dir / "counts/remaining_counts.csv"
    )
    dev = _table_4a("sample_level/sample_qc.csv")
    try:
        assert_frame_equal(legacy.reindex(dev.index), dev, check_dtype=False)
        typer.secho("Table 4a Matches Exactly", fg=GREEN)
    except AssertionError:
        typer.secho("Please carefully check this table\n", fg=YELLOW)
        typer.secho("\n## Legacy\n", fg=MAGENTA)
        typer.secho(legacy.reindex(dev.index).to_string())
        typer.secho("\n## DEV\n", fg=MAGENTA)
        typer.secho(dev.to_string())


def _legacy_table_4b(exclusions):
    return (
        pd.read_csv(exclusions)
        .set_index("Reason")
        .rename_axis("Filter Reason/Description")
        .reindex(["SexDiscordant", "UnexpectedReplicates", "AutosomalHet"])
        .rename(
            index={
                "SexDiscordant": "Sex Discordant",
                "UnexpectedReplicates": "Unexpected Replicates",
                "AutosomalHet": "Autosomal Het",
            },
            columns={"Total": "All Subjects"},
        )
    )


def _table_4b(subject_qc_csv, population_qc_csv):
    from cgr_gwas_qc.reporting.qc_exclusions import subject_exclusion_counts

    return subject_exclusion_counts(subject_qc_csv, population_qc_csv)


def compare_4b(legacy_dir: Path):
    legacy = _legacy_table_4b(legacy_dir / "counts/exclusion_counts.csv")
    dev = _table_4b("subject_level/subject_qc.csv", "subject_level/population_qc.csv")
    try:
        assert_frame_equal(legacy.reindex(dev.index), dev, check_dtype=False)
        typer.secho("Table 4b Matches Exactly", fg=GREEN)
    except AssertionError:
        typer.secho("Please carefully check this table\n", fg=YELLOW)
        typer.secho("\n## Legacy\n", fg=MAGENTA)
        typer.secho(legacy.reindex(dev.index).to_string())
        typer.secho("\n## DEV\n", fg=MAGENTA)
        typer.secho(dev.to_string())


if __name__ == "__main__":
    app()
