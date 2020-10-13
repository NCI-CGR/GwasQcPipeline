from pathlib import Path

import typer

import cgr_gwas_qc
import cgr_gwas_qc.yaml
from cgr_gwas_qc.models import schema_to_dict
from cgr_gwas_qc.models.config import Config

app = typer.Typer(add_completion=False)


@app.command()
def main(
    project_name: str = typer.Option(..., prompt="Project Name", help="The current project title."),
    sample_sheet: Path = typer.Option(
        ..., prompt="Path to LIMs Sample Sheet", exists=True, readable=True, file_okay=True
    ),
):
    """Creates a Gwas Qc Pipeline config file in the current working directory."""
    cfg = schema_to_dict(Config.schema())
    cfg["project_name"] = project_name
    cfg["pipeline_version"] = cgr_gwas_qc.__version__
    cfg["sample_sheet"] = sample_sheet.expanduser().absolute().as_posix()

    cgr_gwas_qc.yaml.write(cfg, "config.yml")


if __name__ == "__main__":
    app()
