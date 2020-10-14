import typer
from snakemake.rules import expand

from cgr_gwas_qc import load_config
from cgr_gwas_qc.validators.gtc import GTCFileError, check_gtc_files
from cgr_gwas_qc.validators.idat import IdatFileError, check_idat_files

app = typer.Typer(add_completion=False)


@app.command()
def main(gtc_check: bool = True, idat_check: bool = True):
    """Pre-flight checks to make sure user input files are readable and complete.

    Pre-flight checks include the Sample Sheet, reference files, Idat files,
    and GTC files. For Idat and GTC files, the user provided file name
    pattern is extended using the columns in the sample sheet.
    """
    cfg = load_config()  # This already validates all Reference files in the config
    typer.echo("Sample Sheet and Reference Files OK.")

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
