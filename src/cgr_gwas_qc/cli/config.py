import logging
import os
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
    sample_sheet: Path = typer.Option(
        ..., prompt="Path to LIMs Sample Sheet", exists=True, readable=True, file_okay=True
    ),
    project_name: Optional[str] = typer.Option(None, help="The current project title."),
    cgems: bool = typer.Option(True, help="Set common default paths for cgems"),
):
    """Creates a Gwas Qc Pipeline config file in the current working directory."""

    ss = SampleSheet(sample_sheet)
    project_name = project_name or ss.header.get("Project Name", "No Project Name").split(";")[0]
    num_samples = ss.data.shape[0]
    snp_array = ss.manifests.get("snp_array", None)
    bpm = ss.manifests.get("bpm", "GSAMD-24v1-0_20011747_A1.bpm")
    output_pattern = get_output_pattern(ss.file_name.stem)

    if cgems:
        cfg = cgems_config(ss.file_name, project_name, num_samples, snp_array, bpm, output_pattern)
    else:
        cfg = general_config(
            ss.file_name, project_name, num_samples, snp_array, bpm, output_pattern
        )

    config_to_yaml(cfg, exclude_none=True)


def cgems_config(
    sample_sheet: Path,
    project_name: str,
    num_samples: int,
    snp_array: Optional[str],
    bpm: str,
    output_pattern: str,
) -> Config:
    bpm_file = Path(f"/DCEG/CGF/Infinium/Resources/Manifests/{bpm}")
    return Config(
        project_name=project_name,
        sample_sheet=sample_sheet,
        num_samples=num_samples,
        snp_array=snp_array,
        num_snps=get_number_snps(bpm_file),
        reference_files=dict(
            illumina_manifest_file=bpm_file,
            thousand_genome_vcf="/DCEG/CGF/Bioinformatics/Production/data/thousG/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz",
            thousand_genome_tbi="/DCEG/CGF/Bioinformatics/Production/data/thousG/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz.tbi",
        ),
        user_files=dict(
            output_pattern=output_pattern,
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


def general_config(
    sample_sheet: Path,
    project_name: str,
    num_samples: int,
    snp_array: Optional[str],
    bpm: str,
    output_pattern: str,
) -> Config:
    bpm_file = Path(f"/path/to/illumina/{bpm}")
    return Config(
        project_name=project_name,
        sample_sheet=sample_sheet,
        num_samples=num_samples,
        snp_array=snp_array,
        num_snps=get_number_snps(bpm_file),
        reference_files=dict(
            illumina_manifest_file=bpm_file,
            thousand_genome_vcf="/path/to/thousand/genome/vcf.gz",
            thousand_genome_tbi="/path/to/thousand/genome/vcf.gz.tbi",
        ),
        user_files=dict(
            output_pattern=output_pattern,
            idat_pattern=dict(
                red="/expample/pattern/wildcards/are/columns/in/sample_sheet/{Project}/{Sample_ID}}_Red.idat",
                green="/expample/pattern/wildcards/are/columns/in/sample_sheet/{Project}/{Sample_ID}}_Grn.idat",
            ),
            gtc_pattern="/expample/pattern/wildcards/are/columns/in/sample_sheet/{Project}/{Sample_ID}.gtc",
        ),
    )


def get_output_pattern(sample_sheet_name):
    if "AnalysisManifest" in sample_sheet_name:
        return "{prefix}/" + sample_sheet_name.replace("AnalysisManifest", "{file_type}").replace(
            ".csv", ".{ext}"
        )

    return "{prefix}/{file_type}.{ext}"


def get_number_snps(manifest_file: Optional[os.PathLike]) -> int:
    if manifest_file:
        try:
            bpm = BeadPoolManifest(manifest_file)
            return bpm.num_loci
        except FileNotFoundError:
            logger.warning(
                "Could not parse the illumina manifest file. Did not set num_snps in the config."
            )
    return 0


if __name__ == "__main__":
    app()
