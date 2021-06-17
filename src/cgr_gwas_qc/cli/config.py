import re
from datetime import date
from pathlib import Path
from textwrap import dedent
from typing import Optional

import pandas as pd
import typer

from cgr_gwas_qc.config import config_to_yaml
from cgr_gwas_qc.models.config import Config, GenomeBuild, Idat, ReferenceFiles, UserFiles
from cgr_gwas_qc.parsers import sample_sheet as sample_sheet_parser
from cgr_gwas_qc.parsers.illumina import BeadPoolManifest

app = typer.Typer(add_completion=False)
TODAY = date.today().strftime("%m%d%Y")

CGEMS_MANIFESTS = Path("/DCEG/CGF/Infinium/Resources/Manifests")
CGEMS_1KG = Path("/DCEG/CGF/Bioinformatics/Production/data/thousG")
CGEMS_SCAN_DATA = Path("/DCEG/CGF/Infinium/ScanData/CGF/ByProject")
CGEMS_QC_DIR = Path("/DCEG/CGF/GWAS/Scans/GSA_Lab_QC")


@app.command()
def main(
    sample_sheet: Path = typer.Option(
        ...,
        "--sample-sheet",
        "-s",
        help="Path to the Sample Sheet",
        prompt="Path to the Sample Sheet",
        exists=True,
        readable=True,
        file_okay=True,
    ),
    project_name: Optional[str] = typer.Option(None, help="The current project title."),
    bpm_file: Optional[Path] = typer.Option(None, help="The BPM File to use."),
    genome_build: GenomeBuild = typer.Option(
        GenomeBuild.hg37, case_sensitive=False, help="Which human genome build to use."
    ),
    include_unused_settings: bool = typer.Option(
        False, "--include-unused-settings", "-u", help="Include unused options in the config."
    ),
    cgems: bool = typer.Option(
        False,
        "--cgems",
        help="Create folder structure and config for a production run on CGEMs/CCAD.",
    ),
    cgems_dev: bool = typer.Option(
        False,
        "--cgems-dev",
        help="Use CGEMs/CCAD standard paths but create the config in the current working directory.",
    ),
    no_prompt: bool = typer.Option(
        False, "-y", help="Don't prompt during CGEMs/CCAD folder creation."
    ),
):
    """Creates a Gwas Qc Pipeline config file in the current working directory."""
    cfg = initialize_config(sample_sheet, project_name, bpm_file, genome_build)

    # Update config to use paths on CGEMs/CCAD or add place holders for other systems
    update_config_for_cgems(cfg) if cgems | cgems_dev else update_config_for_general(cfg)

    cfg = Config(**cfg.dict())

    # Decide where to save the config
    config_file = cgems_production_config_file(cfg, no_prompt) if cgems else "config.yml"

    # Save config
    config_to_yaml(cfg, yaml_file=config_file, exclude_none=not include_unused_settings)
    typer.secho(
        dedent(
            f"""
            Successfully created:
              - `{config_file}`

            Please check and make any necessary changes. Before running the
            workflow, you must run `cgr pre-flight` to check the reference and
            user files and generate the `cgr_sample_sheet.csv` file."
            """
        ),
        fg=typer.colors.GREEN,
    )


def initialize_config(
    sample_sheet: Path,
    project_name: Optional[str],
    bpm_file: Optional[Path],
    genome_build: GenomeBuild,
) -> Config:
    """Initialize config object.

    Here I am using the construct methods on the pydantic data models to skip
    validation. This allows me to incrementally build the config object with
    place holders that are not valid values.
    """
    return Config.construct(
        project_name=project_name,
        sample_sheet=sample_sheet.resolve(),
        genome_build=genome_build,
        num_snps=0,
        reference_files=ReferenceFiles.construct(
            illumina_manifest_file=bpm_file,
            thousand_genome_vcf="/path/to/thousand/genome/vcf.gz",
            thousand_genome_tbi="/path/to/thousand/genome/vcf.gz.tbi",
        ),
        user_files=UserFiles(
            output_pattern="{prefix}/{file_type}.{ext}",
            idat_pattern=Idat(
                red="/expample/pattern/wildcards/are/columns/in/sample_sheet_file/{Project}/{Sample_ID}_Red.idat",
                green="/expample/pattern/wildcards/are/columns/in/sample_sheet_file/{Project}/{Sample_ID}_Grn.idat",
            ),
            gtc_pattern="/expample/pattern/wildcards/are/columns/in/sample_sheet_file/{Project}/{Sample_ID}.gtc",
        ),
    )


def update_config_for_cgems(cfg: Config):
    _add_project_name(cfg)
    _add_sample_sheet_metadata(cfg)
    _add_cgems_reference_files(cfg)
    _add_cgems_user_file_patterns(cfg)
    _update_number_of_snps(cfg)


def update_config_for_general(cfg: Config):
    _add_project_name(cfg)
    _add_sample_sheet_metadata(cfg)
    _update_number_of_snps(cfg)


def cgems_production_config_file(cfg: Config, no_prompt: bool):
    run_dir = _get_cgems_production_run_dir(cfg.project_name)
    if _create_cgems_production_run_dir(run_dir, no_prompt):
        return run_dir / "config.yml"

    typer.secho(
        dedent(
            f"""
            Warning: We are not using:
              - {run_dir}
            Instead, we created the `config.yml` in the current working directory.
            """
        ),
        fg=typer.colors.YELLOW,
    )
    return "config.yml"


def _add_project_name(cfg: Config):
    if cfg.project_name is None:
        cfg.project_name = _sample_sheet_name_to_project_name(cfg.sample_sheet.stem)


def _sample_sheet_name_to_project_name(stem: str):
    """Convert a CGEMs sample sheet name to the project name."""
    if (m := re.match(r"(?P<prefix>.*)_AnalysisManifest_(?P<suffix>.*)", stem)) :
        return "_".join(m.groups())
    return f"GwasQcPipeline_{TODAY}"


def _add_sample_sheet_metadata(cfg: Config):
    if sample_sheet_parser.is_sample_manifest(cfg.sample_sheet):
        ss = sample_sheet_parser.SampleManifest(cfg.sample_sheet)
        cfg.num_samples = ss.data.shape[0]
        cfg.snp_array = ss.manifests.get("snp_array", None)

        if not cfg.reference_files.illumina_manifest_file:
            bpm = ss.manifests.get("bpm", "GSAMD-24v1-0_20011747_A1.bpm")
            cfg.reference_files.illumina_manifest_file = CGEMS_MANIFESTS / bpm
    else:
        df = pd.read_csv(cfg.sample_sheet)
        cfg.num_samples = df.shape[0]


def _add_cgems_reference_files(cfg: Config):
    if cfg.genome_build == "hg37":
        cfg.reference_files.thousand_genome_vcf = (
            CGEMS_1KG / "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz"
        )
        cfg.reference_files.thousand_genome_tbi = (
            CGEMS_1KG / "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz.tbi"
        )
    elif cfg.genome_build == "hg38":
        cfg.reference_files.thousand_genome_vcf = (
            CGEMS_1KG / "hg38/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz"
        )
        cfg.reference_files.thousand_genome_tbi = (
            CGEMS_1KG / "hg38/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz.tbi"
        )
    else:
        typer.secho("We only support genome builds [hg37, hg38]", fg=typer.colors.RED)
        raise typer.Exit(1)


def _add_cgems_user_file_patterns(cfg: Config):
    if "AnalysisManifest" in cfg.sample_sheet.stem:
        file_type = cfg.sample_sheet.stem.replace("AnalysisManifest", "{file_type}")
        cfg.user_files.output_pattern = f"{{prefix}}/{file_type}.{{ext}}"

    cfg.user_files.idat_pattern = Idat(
        red=(
            CGEMS_SCAN_DATA
            / "{Project}/{SentrixBarcode_A}/{SentrixBarcode_A}_{SentrixPosition_A}_Red.idat"
        ).as_posix(),
        green=(
            CGEMS_SCAN_DATA
            / "{Project}/{SentrixBarcode_A}/{SentrixBarcode_A}_{SentrixPosition_A}_Grn.idat"
        ).as_posix(),
    )

    cfg.user_files.gtc_pattern = (
        CGEMS_SCAN_DATA / "{Project}/{SentrixBarcode_A}/{SentrixBarcode_A}_{SentrixPosition_A}.gtc"
    ).as_posix()


def _update_number_of_snps(cfg: Config):
    if cfg.reference_files.illumina_manifest_file.exists():
        try:
            cfg.num_snps = BeadPoolManifest(cfg.reference_files.illumina_manifest_file).num_loci
        except Exception:
            typer.secho(
                "We could not parse the number of SNPs from theillumina manifest file:\n"
                "  - {cfg.reference_files.illumina_manifest_file}",
                fg=typer.colors.YELLOW,
            )


def _get_cgems_production_run_dir(project_name: str) -> Path:
    project_dir = CGEMS_QC_DIR / f"{project_name}/builds"
    num_previous_runs = len(list(project_dir.glob("*")))
    run_number = num_previous_runs + 1
    return project_dir / f"QC_v{run_number}_{TODAY}"


def _create_cgems_production_run_dir(run_dir: Path, no_prompt: bool = False) -> bool:
    create_dir = True if no_prompt else typer.confirm(f"Do you want to create {run_dir}?")

    if create_dir:
        try:
            run_dir.mkdir(parents=True)
            typer.secho(f"Created {run_dir}", fg=typer.colors.GREEN)
            return True
        except FileExistsError:
            typer.secho(f"The directory {run_dir} already exists.", fg=typer.colors.RED)

    return False


if __name__ == "__main__":
    app()
