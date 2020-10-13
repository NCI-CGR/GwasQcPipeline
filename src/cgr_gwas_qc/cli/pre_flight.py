import typer
from snakemake.rules import expand

from cgr_gwas_qc import load_config
from cgr_gwas_qc.data_validation.gtc import GTCFileError, check_gtc_files
from cgr_gwas_qc.data_validation.idat import IdatFileError, check_idat_files

app = typer.Typer(add_completion=False)


@app.command()
def main(gtc_check: bool = True, idat_check: bool = True):
    """Run pre-flight checks to make sure all inputs are safe."""
    cfg = load_config()  # This already does pre-flight checks on Reference files in the config
    typer.echo("Reference Files OK.")

    if idat_check:
        try:
            check_idat_files(
                expand(cfg.config.user_data_patterns.idat.green, **cfg.ss.as_dict("list"))
            )
            check_idat_files(
                expand(cfg.config.user_data_patterns.idat.red, **cfg.ss.as_dict("list"))
            )
            typer.echo("Idat Files OK.")
        except IdatFileError:
            typer.Exit(1)

    if gtc_check:
        try:
            check_gtc_files(expand(cfg.config.user_data_patterns.gtc, **cfg.ss.as_dict("list")))
            typer.echo("GTC Files OK.")
        except GTCFileError:
            typer.Exit(1)


if __name__ == "__main__":
    app()
