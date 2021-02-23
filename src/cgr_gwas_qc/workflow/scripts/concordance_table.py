from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False)


@app.command()
def main(filename: Path, outname: Path):
    (
        pd.read_csv(filename, delim_whitespace=True)
        .assign(concordance=lambda x: x.IBS2 / (x.IBS0 + x.IBS1 + x.IBS2))
        .reindex(["IID1", "IID2", "PI_HAT", "concordance"], axis=1)
        .to_csv(outname)
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        main(Path(snakemake.input[0]), Path(snakemake.output[0]))  # type: ignore # noqa
    else:
        app()
