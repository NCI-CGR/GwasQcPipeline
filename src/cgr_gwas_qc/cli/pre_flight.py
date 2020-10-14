from logging import getLogger
from pathlib import Path

import typer
from snakemake.rules import expand

from cgr_gwas_qc import load_config
from cgr_gwas_qc.validators.gtc import validate as gtc_validate

app = typer.Typer(add_completion=False)
logger = getLogger(__name__)


@app.command()
def main(gtc_check: bool = True, idat_check: bool = False):
    """Pre-flight checks to make sure user input files are readable and complete.

    Pre-flight checks include the Sample Sheet, reference files, Idat files,
    and GTC files. For Idat and GTC files, the user provided file name
    pattern is extended using the columns in the sample sheet.
    """
    cfg = load_config()  # This already validates all Reference files in the config
    typer.echo("Sample Sheet and Reference Files OK.")

    if idat_check:
        logger.warn("Idat pre-flight is not impelmented yet.")

    if gtc_check:
        for gtc in expand(cfg.config.user_data_patterns.gtc, **cfg.ss.to_dict("list")):
            gtc_validate(Path(gtc))
        typer.echo("GTC Files OK.")


if __name__ == "__main__":
    app()
