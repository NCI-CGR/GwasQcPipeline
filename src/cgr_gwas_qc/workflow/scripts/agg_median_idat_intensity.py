from pathlib import Path
from typing import List

import pandas as pd
import typer

app = typer.Typer(add_completion=False)


@app.command()
def main(inputs: List[Path], outfile: Path):
    pd.concat([pd.read_csv(file_name) for file_name in inputs]).to_csv(outfile, index=False)


if __name__ == "__main__":
    if "snakemake" in locals():
        main(
            inputs=[Path(f) for f in snakemake.input],  # type: ignore # noqa
            outfile=Path(snakemake.output[0]),  # type: ignore # noqa
        )
    else:
        app()
