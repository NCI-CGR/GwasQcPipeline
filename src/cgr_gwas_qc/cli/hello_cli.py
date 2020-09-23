"""This is just a test CLI to make sure things are working."""
import typer

app = typer.Typer()


@app.command("hello_world")
def print_hello(name: str = "CGR Bioinformatician"):
    typer.echo(f"Hello there {name}.")


def main():
    app()
