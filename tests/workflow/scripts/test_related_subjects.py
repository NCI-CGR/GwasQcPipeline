from textwrap import dedent

from cgr_gwas_qc.workflow.scripts.related_subjects import (
    create_prune_list,
    merge_overlapping_subjects,
)


def test_merge_overlapping_subjects():
    # GIVEN: a list of subject ID tuples
    subjects = [
        ("S001", "S002"),
        ("S001", "S003"),
        ("S002", "S004"),
        ("S005", "S006"),
    ]

    # WHEN: I combine tuples that have an overlapping ID
    observed = merge_overlapping_subjects(subjects)

    # THEN: I should get a nested list where overlapping IDs are grouped
    # together.
    expected = [["S001", "S002", "S003", "S004"], ["S005", "S006"]]
    assert expected == observed


def test_create_prune_list(tmp_path):
    # GIVEN: A fake plink IMISS file with only the required columns.
    imiss = tmp_path / "subjects.imiss"
    imiss.write_text(
        dedent(
            """\
            IID  F_MISS
            S001  0.01
            S002  0.1
            S003  0.1
            S004  0.1
            S005  0.1
            S006  0.01
            """
        )
    )  # NOTE: call_rate is `1 - F_MISS`

    # and a nested list of grouped subject IDs
    grouped_subjects = [["S001", "S002", "S003", "S004"], ["S005", "S006"]]

    # WHEN: I remove the "best" sample (based on call rate) from each group and
    # create a pruning list.
    observed = create_prune_list(grouped_subjects, imiss)

    # THEN: I should prune these 4 samples.
    expected = ["S002 S002", "S003 S003", "S004 S004", "S005 S005"]
    assert expected == observed
