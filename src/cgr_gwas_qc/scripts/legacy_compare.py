#!/usr/bin/env python
"""Helper script to compare results with the legacy workflow. """
from collections import namedtuple
from functools import partial
from pathlib import Path
from typing import Any, Set

import pandas as pd
import typer
from pandas.testing import assert_series_equal
from scipy.stats import pearsonr
from typer.colors import BLUE, BRIGHT_CYAN, GREEN, MAGENTA, RED, WHITE, YELLOW

from cgr_gwas_qc import load_config, yaml
from cgr_gwas_qc.models.config import Config
from cgr_gwas_qc.testing import chdir, comparison

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
    mix_atol: float = typer.Argument(
        0.001, help="The absolute tolerance level when comparing verifyIDintensities %Mix column."
    ),
    llk_atol: float = typer.Argument(
        500.0, help="The absolute tolerance level when comparing verifyIDintensities LLK column."
    ),
    llk0_atol: float = typer.Argument(
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
            typer.Exit(0)

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
    for cmp in files:
        try:
            cmp.func(cmp.legacy, cmp.current)
            typer.secho(f"{cmp.message} Files Match", fg=GREEN)
        except AssertionError:
            typer.secho(f"{cmp.message} did not match ({cmp.legacy} vs {cmp.current}).", fg=RED)
        except NotImplementedError:
            typer.secho(
                f"Cannot currently compare {cmp.message}",
                fg=YELLOW,
            )


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
    for cmp in files:
        try:
            cmp.func(cmp.legacy, cmp.current)
            typer.secho(f"{cmp.message} Files Match", fg=GREEN)
        except AssertionError:
            typer.secho(f"{cmp.message} did not match ({cmp.legacy} vs {cmp.current}).", fg=RED)
        except NotImplementedError:
            typer.secho(
                f"Cannot currently compare {cmp.message}",
                fg=YELLOW,
            )


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
    for cmp in files:
        try:
            cmp.func(cmp.legacy, cmp.current)
            typer.secho(f"{cmp.message} Files Match", fg=GREEN)
        except AssertionError:
            typer.secho(f"{cmp.message} did not match ({cmp.legacy} vs {cmp.current}).", fg=RED)
        except NotImplementedError:
            typer.secho(
                f"Cannot currently compare {cmp.message}",
                fg=YELLOW,
            )


def compare_concordance(legacy_dir: Path, maf: float, ld: float):
    typer.secho("\n## Comparing Sample Concordance", fg=MAGENTA)

    files = [
        Comparison(
            legacy_dir / "ibd/samples.genome",
            f"sample_level/call_rate_2/samples_maf{maf}_ld{ld}.genome",
            "Sample Level IBD",
            comparison.file_hashes_equal,
        ),
    ]
    for cmp in files:
        try:
            cmp.func(cmp.legacy, cmp.current)
            typer.secho(f"{cmp.message} Files Match", fg=GREEN)
        except AssertionError:
            typer.secho(f"{cmp.message} did not match ({cmp.legacy} vs {cmp.current}).", fg=RED)


def compare_contamination(legacy_dir: Path, mix_atol: float, llk_atol: float, llk0_atol: float):
    typer.secho("\n## Comparing Contamination (Fuzzy Match)", fg=MAGENTA)
    files = [
        Comparison(
            legacy_dir / "all_contam/contam.csv",
            "sample_level/contamination/verifyIDintensity.csv",
            "verifyIDintensity Results",
            partial(
                comparison.assert_verifyIDintensity_equal,
                mix_atol=mix_atol,
                llk_atol=llk_atol,
                llk0_atol=llk0_atol,
            ),
        ),
    ]
    for cmp in files:
        try:
            cmp.func(cmp.legacy, cmp.current)
            typer.secho(f"{cmp.message} Files Fuzzy Match", fg=GREEN)
        except AssertionError:
            typer.secho(
                f"{cmp.message} did not match at current tolerance rates ({cmp.legacy} vs {cmp.current}).",
                fg=RED,
            )

    typer.secho(
        (
            "Fuzzy match thresholds:\n"
            f"  - %Mix +/- {mix_atol}\n"
            f"  - LLK +/- {llk_atol}\n"
            f"  - LLK0 +/- {llk0_atol}"
        ),
        fg=YELLOW,
    )


def _parse_snpweights(filename):
    return (
        pd.read_csv(filename)
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
        pd.read_csv(filename, sep="\t")
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
        table = pd.concat([legacy_counts, dev_counts], axis=1).assign(
            Difference=lambda x: (x.iloc[:, 0] - x.iloc[:, 1]).abs()
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
        pd.read_csv(filename)
        .set_index("Sample_ID")
        .sort_index()
        .replace({"N": False, "Y": True})
        .assign(is_missing_idats=lambda x: ~x.IdatsInProjectDir.astype(bool))
        .astype({**{k: bool for k in exclusion_criteria.keys()}})
        .assign(Call_Rate_2_filter=lambda x: x.Call_Rate_2_filter ^ x.Call_Rate_1_filter)
        .assign(Legacy_Exclusions=lambda x: _get_reason(x, exclusion_criteria))
        .Legacy_Exclusions
    )


def _sample_level_exclusions(filename):
    return (
        pd.read_csv(filename)
        .set_index("Sample_ID")
        .sort_index()
        .analytic_exclusion_reason.rename("Production_Exclusions")
    )


def compare_sample_analytic_exclusions(legacy_dir: Path):
    typer.secho("\n## Sample Exclusion Summary", fg=MAGENTA)
    legacy = _legacy_sample_level_exclusions(legacy_dir / "all_sample_qc.csv")
    dev = _sample_level_exclusions("sample_level/sample_qc.csv")
    try:
        assert_series_equal(legacy, dev, check_names=False)
        typer.secho("Sample Analytic Exclusions Match", fg=GREEN)
    except AssertionError:
        typer.secho("Please carefully check this table\n", fg=YELLOW)
        table = pd.crosstab(legacy, dev, margins=True, margins_name="Total").to_string()
        typer.secho(table)


def _legacy_selected_subjects(filename):
    return pd.read_csv(filename).set_index("Subject_ID").sort_index().squeeze().rename("Legacy")


def _selected_subjects(filename):
    df = (
        pd.read_csv(filename)
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
        assert_series_equal(legacy, dev, check_names=False)
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
        pd.read_csv(filename)
        .set_index("Sample_ID")
        .sort_index()
        .replace({"N": False, "Y": True})
        .astype({**{k: bool for k in exclusion_criteria.keys()}})
        .assign(Legacy_Exclusions=lambda x: _get_reason(x, exclusion_criteria))
        .Legacy_Exclusions
    )


def _subject_level_exclusions(filename):
    from cgr_gwas_qc.workflow.scripts.sample_qc_table import _get_reason

    exclusion_criteria = {
        "is_sex_discordant": "Sex Discordance",
        "is_unexpected_replicate": "Unexpected Replicate",
    }

    return (
        pd.read_csv(filename)
        .set_index("Sample_ID")
        .sort_index()
        .astype({**{k: bool for k in exclusion_criteria.keys()}})
        .assign(Production_Exclusions=lambda x: _get_reason(x, exclusion_criteria))
        .Production_Exclusions
    )


def compare_subject_analytic_exclusions(legacy_dir: Path):
    typer.secho("\n## Subject Exclusion Summary", fg=MAGENTA)
    legacy = _legacy_subject_level_exclusions(legacy_dir / "all_sample_qc.csv")
    dev = _subject_level_exclusions("subject_level/subject_qc.csv")
    try:
        assert_series_equal(legacy.reindex(dev.index), dev, check_names=False)
        typer.secho("Subject Analytic Exclusions Match", fg=GREEN)
    except AssertionError:
        typer.secho("Please carefully check this table\n", fg=YELLOW)
        table = pd.crosstab(
            legacy.reindex(dev), dev, margins=True, margins_name="Total"
        ).to_string()
        typer.secho(table)


if __name__ == "__main__":
    app()
