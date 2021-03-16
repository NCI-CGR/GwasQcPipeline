from pathlib import Path

import pandas as pd


def read_het(filename: Path) -> pd.DataFrame:
    """Parse PLINK's het file format.

    Returns:
        pd.DataFrame:
            A (n x 5) table with the following columns

            .. csv-table::
                :header: name, dtype, description

                **ID** (*index*), object, Sample or Subject ID
                O(HOM), int, Observed number of homozygotes
                E(HOM), int, Expected number of homozygotes
                N(NM), int, Number of non-missing autosomal genotypes
                F, float, Method-of-moments F coefficient estimate

    References:
        - https://www.cog-genomics.org/plink/1.9/basic_stats#ibc
        - https://www.cog-genomics.org/plink/1.9/formats#het


    """
    return (
        pd.read_csv(filename, delim_whitespace=True)
        .drop("FID", axis=1)
        .rename({"IID": "ID"}, axis=1)
        .set_index("ID")
    )
