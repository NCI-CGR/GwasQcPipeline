from datetime import datetime
from pathlib import Path
from typing import Dict, Text

import typer

from cgr_gwas_qc.config import Config
from cgr_gwas_qc.parsers.sample_sheet import SampleSheet
from cgr_gwas_qc.reporting import ExclusionTables, SampleQC, SubjectQC, env
from cgr_gwas_qc.workflow.scripts.known_concordant_samples import (
    read_known_concordance_table,
    read_unknown_concordance_table,
)
from cgr_gwas_qc.workflow.scripts.population_qc_table import read_population_qc
from cgr_gwas_qc.workflow.scripts.sample_qc_table import read_sample_qc
from cgr_gwas_qc.workflow.scripts.snp_qc_table import read_snp_qc

app = typer.Typer(add_completion=False)


@app.command()
def main(
    config: Config,
    sample_sheet_csv: Path,
    snp_qc_csv: Path,
    sample_qc_csv: Path,
    population_qc_csv: Path,
    control_replicates_csv: Path,
    study_replicates_csv: Path,
    unexpected_replicates_csv: Path,
    call_rate_png: Path,
    chrx_inbreeding_png: Path,
    ancestry_png: Path,
    autosomal_heterozygosity_png_dir: Path,
    pca_png_dir: Path,
    hwe_png_dir: Path,
    outfile: Path,
):
    sample_sheet = (
        SampleSheet(sample_sheet_csv)
        .add_group_by_column(config.workflow_params.subject_id_to_use)
        .data
    )
    snp_qc = read_snp_qc(snp_qc_csv)
    sample_qc = read_sample_qc(sample_qc_csv)
    population_qc = read_population_qc(population_qc_csv)
    control_replicates = read_known_concordance_table(control_replicates_csv)
    study_replicates = read_known_concordance_table(study_replicates_csv)
    unexpected_replicates = read_unknown_concordance_table(unexpected_replicates_csv)

    payload = {
        "excel_file_name": Path(
            config.user_files.output_pattern.format(prefix="", file_type="QC_Report", ext="xlsx")
        ).name,  # SR0000-001_1_QC_Report_999999999999.xls
        "date": datetime.now().strftime(
            "%a %b %d %H:%M:%S %Y"
        ),  # https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes
        "config": config,
        "sample_qc": SampleQC.construct(
            config,
            sample_sheet,
            snp_qc,
            sample_qc,
            control_replicates,
            study_replicates,
            unexpected_replicates,
            call_rate_png,
            chrx_inbreeding_png,
            ancestry_png,
        ),
        "subject_qc": SubjectQC.construct(
            population_qc, autosomal_heterozygosity_png_dir, pca_png_dir, hwe_png_dir
        ),
        "exclusion_tables": ExclusionTables.construct(
            config, sample_sheet, sample_qc, population_qc
        ),
    }

    report = create_report(payload)
    outfile.write_text(report)


def create_report(payload: Dict) -> Text:
    template = env.get_template("qc_report.md.j2")
    return template.render(**payload)


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({k: v for k, v in snakemake.params.items()})  # type: ignore # noqa
        defaults["outfile"] = Path(snakemake.output[0])  # type: ignore # noqa
        main(**defaults)
    else:
        app()
