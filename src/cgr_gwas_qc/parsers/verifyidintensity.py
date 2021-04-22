import re
from io import StringIO
from pathlib import Path

import pandas as pd

from cgr_gwas_qc.typing import PathLike


def read(filename: PathLike) -> pd.DataFrame:
    """Parse verifyIDintensity contamination output file.

    The summary table provided by verifyIDintensity is poorly formatted. It
    has comments at the top of the file and uses a row of ``---`` to
    delineate the header row. We parse this file into a DataFrame removing
    comments and rename the ``ID`` column to the ``Sample_ID``. We pull the
    current ``{Sample_ID}`` out of the file name and set it in the table for
    easy aggregation.

    Args:
        filename (Path): Path to a `*/{Sample_ID}.contam.out` file.

    Returns:
        pd.DataFrame:

        - Sample_ID
        - %Mix
        - LLK
        - LLK0
    """
    filepath = Path(filename)
    data = _remove_header_delim(filepath.read_text())

    # Here we use ``delim_whitespace`` because these are multiple spaces and
    # not a tab.
    contam = pd.read_csv(StringIO(data), skiprows=2, delim_whitespace=True)

    # Add the Sample_ID so we can aggregate all samples together.
    sample_id = filepath.stem.split(".")[0]
    contam.ID = sample_id

    return contam.rename({"ID": "Sample_ID"}, axis=1)


def _remove_header_delim(data: str) -> str:
    """Remove the table header delimiter from verifyIDintensity output.

    verifyIDintensity adds a row of "---" to deliminate the header. This
    makes it harder to read it into a DataFrame.

    Returns:
        str: The data without the "---" row.
    """
    return re.sub(r"\n-{2,}\n", "\n", data)
