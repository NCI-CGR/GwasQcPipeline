import shutil
from datetime import date
from pathlib import Path
from textwrap import dedent

import typer

from cgr_gwas_qc.cli.config import (
    CGEMS_QC_DIR,
    _create_cgems_production_run_dir,
    _get_cgems_production_run_dir,
)

app = typer.Typer(add_completion=False)
TODAY = date.today().strftime("%m%d%Y")


@app.command()
def main(
    project_dir: Path = typer.Argument(
        ..., help="Path to the CGEMs/CCAD project you want to create a new run for."
    )
):
    """Copy settings to a new CGEMs/CCAD production QC run.


    This command is only meant to be used on CGEMs/CCAD. It copies the
    config.yml and cgr_sample_sheet.csv to a new run directory.
    """
    project_name = _get_project_name(project_dir)
    previous_run_dir = _get_previous_cgems_production_run_dir(project_name)
    run_dir = _get_cgems_production_run_dir(project_name)
    if not _create_cgems_production_run_dir(run_dir, no_prompt=True):
        typer.secho(
            f"Failed to create {run_dir}. Please check paths and try again.", fg=typer.colors.RED
        )
        raise typer.Exit(1)

    _clone_previous_run(previous_run_dir, run_dir)

    typer.secho(
        dedent(
            f"""
            Successfully
              - created {run_dir}
              - copied config.yml
              - copied cgr_sample_sheet.csv

            {run_dir} is ready to be submitted
            """
        ),
        fg=typer.colors.GREEN,
    )


def _get_project_name(project_dir: Path) -> str:
    """Figure out what path the user provided and return the project name

    The user could provide a number of possible paths:
      - project level: /DCEG/CGF/GWAS/Scans/GSA_Lab_QC/<project>
      - builds level: /DCEG/CGF/GWAS/Scans/GSA_Lab_QC/<project>/builds
      - run level: /DCEG/CGF/GWAS/Scans/GSA_Lab_QC/<project>/builds/<run>

    Returns:
        str: <project> from above
    """
    if (project_dir / "builds").exists():
        # project level given
        return project_dir.name

    if project_dir.name == "builds":
        # builds level given
        return project_dir.parent.name

    if project_dir.parent.name == "builds":
        # run level given
        return project_dir.parents[1].name

    typer.secho("Could not find the `builds` directory.", fg=typer.colors.RED)
    raise typer.Exit(1)


def _get_previous_cgems_production_run_dir(project_name: str) -> Path:
    builds_dir = CGEMS_QC_DIR / f"{project_name}/builds"
    return sorted(builds_dir.glob("*"))[-1]


def _clone_previous_run(previous_run_dir: Path, run_dir: Path):
    shutil.copyfile(previous_run_dir / "config.yml", run_dir / "config.yml")
    shutil.copyfile(previous_run_dir / "cgr_sample_sheet.csv", run_dir / "cgr_sample_sheet.csv")


if __name__ == "__main__":
    app()
