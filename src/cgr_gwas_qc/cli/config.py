import os
from pathlib import Path
from typing import Optional

import typer

from cgr_gwas_qc.config import config_to_yaml
from cgr_gwas_qc.models.config import Config, GenomeBuild
from cgr_gwas_qc.parsers.illumina import BeadPoolManifest
from cgr_gwas_qc.parsers.sample_sheet import SampleManifest
from cgr_gwas_qc.typing import PathLike

app = typer.Typer(add_completion=False)


@app.command()
def main(
    sample_sheet: Path = typer.Option(
        ..., prompt="Path to LIMs Sample Sheet", exists=True, readable=True, file_okay=True
    ),
    project_name: Optional[str] = typer.Option(None, help="The current project title."),
    bpm_file: Optional[Path] = typer.Option(None, help="The BPM File to use."),
    cgems: bool = typer.Option(True, help="Set common default paths for cgems"),
    genome_build: GenomeBuild = typer.Option(
        GenomeBuild.hg37, case_sensitive=False, help="Which human genome build to use."
    ),
):
    """Creates a Gwas Qc Pipeline config file in the current working directory."""

    ss = SampleManifest(sample_sheet)
    project_name = project_name or ss.header.get("Project Name", "No Project Name").split(";")[0]
    num_samples = ss.data.shape[0]
    snp_array = ss.manifests.get("snp_array", None)
    output_pattern = get_output_pattern(ss.file_name.stem)

    bpm = (
        bpm_file.resolve()
        if bpm_file is not None
        else ss.manifests.get("bpm", "GSAMD-24v1-0_20011747_A1.bpm")
    )

    if cgems:
        cfg = cgems_config(
            ss.file_name,
            project_name,
            genome_build.value,
            num_samples,
            snp_array,
            bpm,
            output_pattern,
        )
    else:
        cfg = general_config(
            ss.file_name,
            project_name,
            genome_build.value,
            num_samples,
            snp_array,
            bpm,
            output_pattern,
        )

    config_to_yaml(cfg, exclude_none=True)
    typer.secho(
        "Please check the generated `config.yml` and make any necessary changes. "
        "Before running the workflow, you must run `cgr pre-flight` to generate "
        "the `cgr_sample_sheet.csv` file.",
        fg=typer.colors.GREEN,
    )


def cgems_config(
    sample_sheet: Path,
    project_name: str,
    genome_build: str,
    num_samples: int,
    snp_array: Optional[str],
    bpm: PathLike,
    output_pattern: str,
) -> Config:
    bpm_file = (
        Path(f"/DCEG/CGF/Infinium/Resources/Manifests/{bpm}") if isinstance(bpm, str) else bpm
    )

    if genome_build == "hg37":
        vcf = "/DCEG/CGF/Bioinformatics/Production/data/thousG/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz"
        tbi = "/DCEG/CGF/Bioinformatics/Production/data/thousG/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz.tbi"
    else:
        vcf = "/DCEG/CGF/Bioinformatics/Production/data/thousG/hg38/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz"
        tbi = "/DCEG/CGF/Bioinformatics/Production/data/thousG/hg38/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz.tbi"

    return Config(
        project_name=project_name,
        sample_sheet=sample_sheet,
        genome_build=genome_build,
        num_samples=num_samples,
        snp_array=snp_array,
        num_snps=get_number_snps(bpm_file),
        reference_files=dict(
            illumina_manifest_file=bpm_file,
            thousand_genome_vcf=vcf,
            thousand_genome_tbi=tbi,
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
    genome_build: str,
    num_samples: int,
    snp_array: Optional[str],
    bpm: PathLike,
    output_pattern: str,
) -> Config:
    bpm_file = Path(f"/path/to/illumina/{bpm}") if isinstance(bpm, str) else bpm
    return Config(
        project_name=project_name,
        sample_sheet=sample_sheet,
        genome_build=genome_build,
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


def get_output_pattern(sample_sheet_stem: str):
    if "AnalysisManifest" in sample_sheet_stem:
        return "{prefix}/" + sample_sheet_stem.replace("AnalysisManifest", "{file_type}") + ".{ext}"

    return "{prefix}/{file_type}.{ext}"


def get_number_snps(manifest_file: Optional[os.PathLike]) -> int:
    if manifest_file:
        try:
            bpm = BeadPoolManifest(manifest_file)
            return bpm.num_loci
        except FileNotFoundError:
            typer.secho(
                f"Note: we could no parse the illumina manifest file:\n  - {manifest_file}",
                fg=typer.colors.YELLOW,
            )
    return 0


if __name__ == "__main__":
    app()
