from io import StringIO
from textwrap import dedent

import pandas as pd
from pandas.testing import assert_frame_equal

from cgr_gwas_qc.parsers.verifyidintensity import _remove_header_delim


def test_remove_header_delim():
    """Remove header delimiter from verifIDintensity output."""
    # GIVEN: a table with the header offset by "----" (like those from verifyIDintensity)
    observed = dedent(
        """\
        ID  %Mix  LLK LKK0
        ------------------
        0   0.95 -2  0.2
        """
    )

    # and our expected out of no "----"
    expected = dedent(
        """\
        ID  %Mix  LLK LKK0
        0   0.95 -2  0.2
        """
    )

    # WHEN-THEN: we parse the table it should remove the row of "----"
    assert _remove_header_delim(observed) == expected


def test_pandas_parsing_multiple_spaces():
    """Test that pandas is parsing strange spacing correctly."""
    # GIVEN: a table separated by a varying number of spaces (like those from verifyIDintensity)
    observed = dedent(
        """\
        ID  %Mix  LLK    LLK0
        0   0.95 -2  0.2
        """
    )

    # and a dataframe version of the data
    expected = pd.DataFrame({"ID": [0], "%Mix": [0.95], "LLK": [-2], "LLK0": [0.2]})

    # WHEN-THEN: running read_csv with `delim_whitespace` should produce the correct table
    assert_frame_equal(pd.read_csv(StringIO(observed), delim_whitespace=True), expected)
