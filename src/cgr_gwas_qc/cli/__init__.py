"""CGR Gwas QC Pipeline Command Line Interface.

Welcome to the CGR Gwas Qc Pipeline CLI. Through this intervace you can try a
demo, create a config file, run pre-flight checks, and interact with
snakemake.
"""
import typer

from cgr_gwas_qc.cli import config, pre_flight, snakemake

app = typer.Typer(add_completion=False, help=__doc__)
app.command("config")(config.main)
app.command("pre-flight")(pre_flight.main)
app.command("snakemake")(snakemake.main)


if __name__ == "__main__":
    app()
