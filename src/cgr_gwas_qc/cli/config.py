from pathlib import Path

import typer

from cgr_gwas_qc.config import config_to_yaml
from cgr_gwas_qc.models import config as cfg_models

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

    config_to_yaml(cfg, exclude_none=True)


def cgems_config(project_name, sample_sheet):
    return cfg_models.Config.construct(
        project_name=project_name,
        sample_sheet=sample_sheet,
        reference_files=cfg_models.ReferenceFiles.construct(
            illumina_manifest_file="/DCEG/CGF/Infinium/Resources/Manifests/GSAMD-Files/build37/GSAMD-24v1-0_20011747_A1.bpm",
            thousand_genome_vcf="/DCEG/CGF/Bioinformatics/Production/data/thousG/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz",
            thousand_genome_tbi="/DCEG/CGF/Bioinformatics/Production/data/thousG/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz.tbi",
        ),
        user_files=cfg_models.UserFiles.construct(
            idat_pattern=cfg_models.Idat.construct(
                red="/DCEG/CGF/Infinium/ScanData/CGF/ByProject/{Project}/{SentrixBarcode_A}/{SentrixBarcode_A}_{SentrixPosition_A}_Red.idat",
                green="/DCEG/CGF/Infinium/ScanData/CGF/ByProject/{Project}/{SentrixBarcode_A}/{SentrixBarcode_A}_{SentrixPosition_A}_Grn.idat",
            ),
            gtc_pattern="/DCEG/CGF/Infinium/ScanData/CGF/ByProject/{Project}/{SentrixBarcode_A}/{SentrixBarcode_A}_{SentrixPosition_A}.gtc",
        ),
        software_params=cfg_models.SoftwareParams(),
        workflow_params=cfg_models.WorkflowParams(),
        env_modules=cfg_models.EnvModules.construct(
            plink2="plink2/1.90b5", eigensoft="eigensoft/7.2.1", r="R/3.4.0",
        ),
    )


def general_config(project_name, sample_sheet):
    return cfg_models.Config.construct(
        project_name=project_name,
        sample_sheet=sample_sheet,
        reference_files=cfg_models.ReferenceFiles.construct(
            illumina_manifest_file="/path/to/illumina/manifest.bpm",
            thousand_genome_vcf="/path/to/thousand/genome/vcf.gz",
            thousand_genome_tbi="/path/to/thousand/genome/vcf.gz.tbi",
        ),
        user_files=cfg_models.UserFiles.construct(
            idat_pattern=cfg_models.Idat.construct(
                red="/expample/pattern/wildcards/are/columns/in/sample_sheet/{Project}/{Sample_ID}}_Red.idat",
                green="/expample/pattern/wildcards/are/columns/in/sample_sheet/{Project}/{Sample_ID}}_Grn.idat",
            ),
            gtc_pattern="/expample/pattern/wildcards/are/columns/in/sample_sheet/{Project}/{Sample_ID}.gtc",
        ),
        software_params=cfg_models.SoftwareParams(),
        workflow_params=cfg_models.WorkflowParams(),
    )


if __name__ == "__main__":
    app()
