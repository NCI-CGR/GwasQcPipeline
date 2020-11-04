from pathlib import Path

import typer

from cgr_gwas_qc.config import create_yaml_config

app = typer.Typer(add_completion=False)


@app.command()
def main(
    project_name: str = typer.Option(..., prompt="Project Name", help="The current project title."),
    sample_sheet: Path = typer.Option(
        ..., prompt="Path to LIMs Sample Sheet", exists=True, readable=True, file_okay=True
    ),
):
    """Creates a Gwas Qc Pipeline config file in the current working directory."""
    create_yaml_config(project_name, sample_sheet.expanduser().absolute())


if __name__ == "__main__":
    app()
