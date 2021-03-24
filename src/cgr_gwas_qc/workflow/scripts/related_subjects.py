#!/usr/bin/env python
from pathlib import Path
from typing import Generator, Optional

import networkx as nx
import pandas as pd
import typer
from numpy.random import RandomState

from cgr_gwas_qc.parsers import plink

app = typer.Typer(add_completion=False)


@app.command()
def main(
    genome: Path = typer.Argument(..., help="Path to the PLINK `.genome` file", exists=True),
    pi_hat_threshold: float = typer.Argument(  # noqa
        ..., help="The threshold to consider subjects related."
    ),
    relatives: Path = typer.Argument(..., help="Path to save a list of subjects to remove."),
    to_remove: Path = typer.Argument(..., help="Path to save a list of subjects to remove."),
):
    """Create a related subjects pruning list.

    The goal of this script is to prune related subjects based on `PI_HAT`.
    It uses a graph based algorithm to select which subjects should be
    pruned. It does this by doing the following:

    1. Uses `PI_HAT` from IBS/IBD estimation to identify pairwise groups of
       related subjects.
    2. Builds a graph where subjects (nodes) are connected by edges if their
       `PI_HAT` was above the pi_hat_threshold.
    3. Identifies subjects that have the most relatives (ie., nodes with max
       degree).
    4. Randomly selects and removes one of these highly connected subjects.
       Adds this subject to the output pruning list.
    5. Removes any subjects that no longer have a relative in the graph
       (i.e., isolated nodes with no edges).
    6. Repeats 3-5 until there are no subjects left in the graph.

    Saves a list of subjects to prune in a format compatible with PLINK's
    `--remove` option.
    """
    pairwise_related_subjects = list(
        plink.read_genome(genome)
        .query("PI_HAT > @pi_hat_threshold")  # Ignore subjects under the pi_hat_threshold
        .reindex(["ID1", "ID2"], axis=1)
        .itertuples(index=False)
    )

    if not pairwise_related_subjects:
        # no related subjects so create empty files
        relatives.write_text("QC_Family_ID,relatives\n")
        to_remove.touch()
        return None

    # Create a graph where subjects are nodes and edges are relatedness at a
    # `PI_HAT > pi_hat_threshold`.
    G = nx.Graph()
    G.add_edges_from(pairwise_related_subjects)

    # Create a list of relatives at the current PI_HAT threshold
    create_qc_families(pairwise_related_subjects).to_csv(relatives)

    # Create a pruning list by removing the most connected subjects first.
    prune_list = create_prune_list(G)
    to_remove.write_text("\n".join(prune_list))


def create_qc_families(G: nx.Graph) -> pd.Series:
    """Create a summary table of putative relatives.

    For each subgraph of relatives, concatenates node IDs into a string and
    assigns an arbitrary `QC_Family_ID`.

    Args:
        G: A graph where nodes are subjects and edges indicate two
          subjects are related.

    Returns:
        pd.Series:

            .. csv-table::
                :header: name, dtype, description

                **QC_Family_ID** (*index*), string, An arbitrary ID assigned to each related subgraph.
                relatives, string, A list of related IDs concatenated together with a `|`.

    """
    return pd.Series(
        {
            f"fam{i}": "|".join(sorted(subgraph))
            for i, subgraph in enumerate(sorted(nx.connected_components(G), key=len), start=1)
        },
        name="relatives",
    ).rename_axis("QC_Family_ID")


def create_prune_list(
    G: nx.Graph, seed: Optional[RandomState] = None
) -> Generator[str, None, None]:
    """Generator creating a list of subjects to prune.

    Keeps running until the graph `G` is empty.

    Args:
        G: A graph where nodes are subjects and edges indicate two
          subjects are related.
        seed: A numpy random seed. Added to allow tests to have consistent
          output.

    Yields:
        Generator: A string in the format "Subject_ID Subject_ID"
    """
    while len(G) > 0:
        yield from _prune(G, seed)


def _prune(G: nx.Graph, seed: Optional[RandomState]) -> Generator[str, None, None]:
    subject_to_remove = (
        pd.Series(dict(G.degree))
        .pipe(lambda x: x[x == x.max()])  # select the subjects with the most relatives
        .sample(n=1, random_state=seed)  # randomly select one of these highly connect subjects
        .index[0]
    )
    G.remove_node(subject_to_remove)  # remove selected subject from the graph

    non_connected_nodes = list(nx.isolates(G))  # Find subjects with no remaining relatives
    G.remove_nodes_from(non_connected_nodes)  # Remove isolated subjects from the graph

    yield f"{subject_to_remove} {subject_to_remove}"  # add subject_to_remove to pruning list


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({"genome": Path(snakemake.input[0])})  # type: ignore # noqa
        defaults.update({k: Path(v) for k, v in snakemake.output.items()})  # type: ignore # noqa
        defaults.update(snakemake.params)  # type: ignore # noqa
        main(**defaults)
    else:
        app()
