from io import StringIO
from textwrap import dedent

import pandas as pd
from pandas.testing import assert_frame_equal


def test_save_subjects_used_for_study(tmp_path):
    from cgr_gwas_qc.workflow.scripts.sample_lists_from_qc_flags import (
        _save_subjects_used_for_study,
    )

    # GIVEN: A slimmed down QC report
    df = pd.read_csv(
        StringIO(
            dedent(
                """
        Group_By_Subject_ID,Sample_ID,Subject_Representative,Subject_Dropped_From_Study,Internal_Control
        SB001,S001,True,False,False
        SB001,S002,False,False,False
        SB002,S003,False,True,False
        SB002,S004,False,True,False
        SB003,I001,False,False,True
        """
            )
        ),
    )

    # WHEN: I save the subjects used for study
    _save_subjects_used_for_study(df, tmp_path / "test.csv")

    # THEN: observed should be the same as expected
    #  - Internal Controls should be removed
    #  - Representative samples should be retained
    #  - Subjects with no representative sample should be Sample_ID = "NA"
    obs_ = pd.read_csv(tmp_path / "test.csv").fillna("NA")
    exp_ = pd.DataFrame({"Subject_ID": ["SB001", "SB002"], "Sample_ID": ["S001", "NA"]})
    assert_frame_equal(obs_, exp_, check_names=False, check_like=True)
