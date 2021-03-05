import logging
from pathlib import Path
from typing import Optional

import typer

from cgr_gwas_qc.config import config_to_yaml
from cgr_gwas_qc.models.config import Config
from cgr_gwas_qc.parsers.illumina import BeadPoolManifest
from cgr_gwas_qc.parsers.sample_sheet import SampleSheet

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
app = typer.Typer(add_completion=False)


@app.command()
def main(
    project_name: str = typer.Option(..., prompt="Project Name", help="The current project title."),
    sample_sheet: Path = typer.Option(
        ..., prompt="Path to LIMs Sample Sheet", exists=True, readable=True, file_okay=True
    ),
    cgems: bool = typer.Option(True, help="Set common default paths for cgems"),
):
    """Creates a Gwas Qc Pipeline config file in the current working directory."""

    if cgems:
        cfg = cgems_config(project_name, sample_sheet)
    else:
        cfg = general_config(project_name, sample_sheet)

    cfg.num_samples = get_number_samples(sample_sheet)
    cfg.num_snps = get_number_snps(cfg.reference_files.illumina_manifest_file)

    config_to_yaml(cfg, exclude_none=True)


def cgems_config(project_name, sample_sheet):
    return Config(
        project_name=project_name,
        sample_sheet=sample_sheet,
        reference_files=dict(
            illumina_manifest_file="/DCEG/CGF/Infinium/Resources/Manifests/GSAMD-Files/build37/GSAMD-24v1-0_20011747_A1.bpm",
            thousand_genome_vcf="/DCEG/CGF/Bioinformatics/Production/data/thousG/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz",
            thousand_genome_tbi="/DCEG/CGF/Bioinformatics/Production/data/thousG/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz.tbi",
        ),
        user_files=dict(
            idat_pattern=dict(
                red="/DCEG/CGF/Infinium/ScanData/CGF/ByProject/{Project}/{SentrixBarcode_A}/{SentrixBarcode_A}_{SentrixPosition_A}_Red.idat",
                green="/DCEG/CGF/Infinium/ScanData/CGF/ByProject/{Project}/{SentrixBarcode_A}/{SentrixBarcode_A}_{SentrixPosition_A}_Grn.idat",
            ),
            gtc_pattern="/DCEG/CGF/Infinium/ScanData/CGF/ByProject/{Project}/{SentrixBarcode_A}/{SentrixBarcode_A}_{SentrixPosition_A}.gtc",
        ),
        env_modules=dict(
            plink2="plink2/1.90b5", eigensoft="eigensoft/7.2.1", r="R/3.4.0", graf="graf/2.3.1"
        ),
    )


def general_config(project_name, sample_sheet):
    return Config(
        project_name=project_name,
        sample_sheet=sample_sheet,
        reference_files=dict(
            illumina_manifest_file="/path/to/illumina/manifest.bpm",
            thousand_genome_vcf="/path/to/thousand/genome/vcf.gz",
            thousand_genome_tbi="/path/to/thousand/genome/vcf.gz.tbi",
        ),
        user_files=dict(
            idat_pattern=dict(
                red="/expample/pattern/wildcards/are/columns/in/sample_sheet/{Project}/{Sample_ID}}_Red.idat",
                green="/expample/pattern/wildcards/are/columns/in/sample_sheet/{Project}/{Sample_ID}}_Grn.idat",
            ),
            gtc_pattern="/expample/pattern/wildcards/are/columns/in/sample_sheet/{Project}/{Sample_ID}.gtc",
        ),
    )


def get_number_samples(sample_sheet: Path) -> Optional[int]:
    try:
        ss = SampleSheet(sample_sheet)
        return ss.data.shape[0]
    except FileNotFoundError:
        logger.warning(
            "Could not parse the sample sheet file. Did not set num_samples in the config."
        )
        return None


def get_number_snps(manifest_file: Path) -> Optional[int]:
    try:
        bpm = BeadPoolManifest(manifest_file)
        return bpm.num_loci
    except FileNotFoundError:
        logger.warning(
            "Could not parse the illumina manifest file. Did not set num_snps in the config."
        )
        return None


if __name__ == "__main__":
    app()
