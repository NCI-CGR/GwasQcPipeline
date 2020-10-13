import typer

from cgr_gwas_qc.cli import config, pre_flight, snakemake

app = typer.Typer(add_completion=False)
app.command("config")(config.main)
app.command("pre-flight")(pre_flight.main)
app.command("snakemake")(snakemake.main)


if __name__ == "__main__":
    app()
