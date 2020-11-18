from io import StringIO
from textwrap import dedent

import pandas as pd
from pandas.testing import assert_frame_equal

from cgr_gwas_qc.parsers.verifyidintensity import remove_header_delim


def test_remove_header_delim():
    """Remove header delimiter from verifIDintensity output."""
    observed = dedent(
        """\
        ID  %Mix  LLK LKK0
        ------------------
        0   0.95 -2  0.2
        """
    )

    expected = dedent(
        """\
        ID  %Mix  LLK LKK0
        0   0.95 -2  0.2
        """
    )

    assert remove_header_delim(observed) == expected


def test_pandas_parsing_multiple_spaces():
    """Test that pandas is parsing strange spacing correctly."""
    observed = dedent(
        """\
        ID  %Mix  LLK    LLK0
        0   0.95 -2  0.2
        """
    )

    expected = pd.DataFrame({"ID": [0], "%Mix": [0.95], "LLK": [-2], "LLK0": [0.2]})

    assert_frame_equal(pd.read_csv(StringIO(observed), delim_whitespace=True), expected)
