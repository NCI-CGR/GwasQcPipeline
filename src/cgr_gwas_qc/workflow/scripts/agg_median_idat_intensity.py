from typing import Sequence

import pandas as pd

from cgr_gwas_qc.typing import PathLike


def main(inputs: Sequence[PathLike], outfile: PathLike):
    if len(inputs) == 1:
        df = pd.read_csv(inputs[0])
    else:
        df = pd.concat([pd.read_csv(filename) for filename in inputs], ignore_index=True)

    df.to_csv(outfile, index=False)


if __name__ == "__main__" and "snakemake" in locals():
    main(
        inputs=list(snakemake.input),  # type: ignore # noqa
        outfile=snakemake.output[0],  # type: ignore # noqa
    )
