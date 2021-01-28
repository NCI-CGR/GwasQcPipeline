#!/usr/bin/env python

from itertools import groupby, product
from pathlib import Path
from typing import Iterable, List, Tuple

import pandas as pd
import typer

app = typer.Typer(add_completion=False)


@app.command()
def main(
    imiss: Path = typer.Argument(..., help="Path to the subject level call rate.", exists=True),
    ibs: Path = typer.Argument(..., help="Path to the PLINK `.genome` file", exists=True),
    pi_hat_threshold: float = typer.Argument(  # noqa
        ..., help="The threshold to consider subjects related."
    ),
    out_file: Path = typer.Argument(..., help="Path to save list subjects to remove."),
):
    """Create related subjects pruning list.

    Identifies groups of subjects that are related based on PI_HAT. For each
    group, keep the subject with the highest call rate. Save the remaining
    subjects to a file for use with PLINK's `--remove` command.
    """
    pairwise_related_subjects = list(
        pd.read_csv(ibs, delim_whitespace=True)
        .query("PI_HAT > @pi_hat_threshold")
        .reindex(["IID1", "IID2"], axis=1)
        .itertuples(index=False)
    )

    if not pairwise_related_subjects:
        # no related subjects so create empty pruning file
        out_file.touch()
        return None

    grouped_related_subjects = merge_overlapping_subjects(pairwise_related_subjects)
    prune_list = create_prune_list(grouped_related_subjects, imiss)
    out_file.write_text("\n".join(prune_list))


def merge_overlapping_subjects(subject_list: Iterable[Tuple[str, str]]) -> List[List[str]]:
    """Create a list of overlapping subjects.

    Given a list of subject ID tuples [(sub1, sub2), (sub1, sub3), (sub4,
    sub5)]. Merge tuples together that have overlapping subject IDs [[sub1,
    sub2, sub3], [sub4, sub5]].
    """
    # Convert Tuples to Sets
    subject_sets = [set(x) for x in subject_list]

    # Do all pairwise comparisons and update each set in place if `a` overlaps
    # with `b`. Update is an in place union.
    for a, b in product(subject_sets, subject_sets):
        if a.intersection(b):
            a.update(b)
            b.update(a)

    # Convert from a list of sets to a list of lists
    combined_lists = sorted([sorted(list(x)) for x in subject_sets])

    # Remove duplicates and return list of lists
    return [grouped_ids for grouped_ids, _ in groupby(combined_lists)]


def create_prune_list(grouped_subjects: List[List[str]], imiss: Path) -> List[str]:
    """Create a list of subjects for pruning.

    For each group of subjects, identify the subject with the highest call
    rate. Then return the remaining subjects to mark for pruning.

    Returns:
        A list of subjects to remove using PLINK's `--remove`.
    """
    subject2call_rate = (
        pd.read_csv(imiss, delim_whitespace=True)
        .set_index("IID")
        .rename_axis("Subject_ID")
        .assign(call_rate=lambda x: 1 - x.F_MISS)
        .call_rate
    )

    subjects_to_prune = []
    for related_subjects in grouped_subjects:
        best_subject_idx = subject2call_rate.reindex(related_subjects).argmax()
        del related_subjects[best_subject_idx]  # don't prune the best subject
        subjects_to_prune.extend(related_subjects)

    return [f"{subID} {subID}" for subID in subjects_to_prune]


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({"out_file": Path(snakemake.output[0])})  # type: ignore # noqa
        defaults.update(snakemake.params)  # type: ignore # noqa
        main(**defaults)
    else:
        app()
