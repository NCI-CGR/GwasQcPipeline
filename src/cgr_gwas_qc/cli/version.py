import typer

from cgr_gwas_qc.version import __version__


def main():
    """Print the current workflow version."""
    typer.echo(__version__)
