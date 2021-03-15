from pathlib import Path

import numpy as np
import pandas as pd


class Eigenvec:
    """Eigensoft eigenvec file parser.

    Attributes:

        filename (Path): Path to the eigenvec file from ``smartpca``.
        components (pd.DataFrame): A (n x 11) table of principal components

            .. csv-table::
                :header: name, dtype, description

                **ID** (*index*), object, Sample or Subject ID
                PC1, object, The first principal componenet (PC1).
                PC2, object, The second principal componenet.
                ..., ..., The third to ninth principal components.
                PC10, object, The tenth principal componenet.

        values (np.ndarray): A vector of eigenvalues for PC1 to PC10.

    References:
        - https://www.hsph.harvard.edu/alkes-price/software/
        - https://www.cell.com/ajhg/fulltext/S0002-9297(16)00003-3
    """

    def __init__(self, filename: Path):
        self.filename = filename
        self.components = self.get_components()
        self.values = self.get_values()

    def get_components(self) -> pd.DataFrame:
        return pd.read_csv(
            self.filename,
            skiprows=1,
            delim_whitespace=True,
            header=None,
            usecols=range(11),
            names=["ID", *(f"PC{i}" for i in range(1, 11))],
        ).set_index("ID")

    def get_values(self) -> np.ndarray:
        return pd.read_csv(
            self.filename,
            nrows=1,
            header=None,
            delim_whitespace=True,
            usecols=range(1, 11),
            names=[f"PC{i}" for i in range(1, 11)],
        ).values[0]
