from pathlib import Path
from typing import Sequence

import pandas as pd

from cgr_gwas_qc.parsers import verifyidintensity
from cgr_gwas_qc.typing import PathLike


def main(inputs: Sequence[PathLike], outfile: PathLike):
    if len(inputs) == 1:
        df = _parse(inputs[0])
    else:
        df = pd.concat([_parse(filename) for filename in inputs], ignore_index=True)

    df.to_csv(outfile, index=False)


def _parse(filename: PathLike):
    path = Path(filename)
    if "---" in path.read_text():
        return verifyidintensity.read(filename)
    return pd.read_csv(filename)


if __name__ == "__main__" and "snakemake" in locals():
    main(
        inputs=list(snakemake.input),  # type: ignore # noqa
        outfile=snakemake.output[0],  # type: ignore # noqa
    )
