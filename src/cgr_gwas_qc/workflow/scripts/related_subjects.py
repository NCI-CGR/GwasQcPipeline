#!/usr/bin/env python

from pathlib import Path
from typing import Generator, Optional

import networkx as nx
import pandas as pd
import typer
from numpy.random import RandomState

app = typer.Typer(add_completion=False)


@app.command()
def main(
    ibs: Path = typer.Argument(..., help="Path to the PLINK `.genome` file", exists=True),
    pi_hat_threshold: float = typer.Argument(  # noqa
        ..., help="The threshold to consider subjects related."
    ),
    out_file: Path = typer.Argument(..., help="Path to save list subjects to remove."),
):
    """Create related subjects pruning list.

    The goal of this script is to prune related subjects based on PI_HAT. It
    uses a greedy based graph algorithm to select which subjects should be
    pruned. It does the following:

    1. Uses PI_HAT from IBS/IBD estimation to identify pairwise groups of
       related subjects.
    2. Builds a graph where subjects (nodes) are connected by edges if their
       PI_HAT was above the pi_hat_threshold.
    3. Identifies subjects that have the most relatives (ie., nodes with max
       degree).
    4. Randomly selects one of these subjects for removal.
    5. Removes any subjects that no longer have a relative in the graph
       (i.e., isolated nodes with no edges).
    6. Repeats 3-5 until there are no subjects left in the graph.

    Writes out a list of subjects to prune in a format compatible with
    PLINK's `--remove` option.
    """
    pairwise_related_subjects = list(
        pd.read_csv(ibs, delim_whitespace=True)
        .query("PI_HAT > @pi_hat_threshold")  # Ignore subjects under the pi_hat_threshold
        .reindex(["IID1", "IID2"], axis=1)
        .itertuples(index=False)
    )

    if not pairwise_related_subjects:
        # no related subjects so create empty pruning file
        out_file.touch()
        return None

    # Create a graph with subjects as nodes.
    G = nx.Graph()
    G.add_edges_from(pairwise_related_subjects)

    # Create pruning list and save.
    prune_list = create_prune_list(G)
    out_file.write_text("\n".join(prune_list))


def create_prune_list(
    G: nx.Graph, seed: Optional[RandomState] = None
) -> Generator[str, None, None]:
    """Generator creating a list of subjects to prune.

    Keeps running until the graph `G` is empty.

    Args:
        G: A graph where nodes are subjects and edges indicate two
          subjects are related.
        seed: A numpy random seed. This is required for testing.

    Yields:
        Generator[str, None, None]: A string in the format "Subject_ID Subject_ID"
    """
    while len(G) > 0:
        yield from _prune(G, seed)


def _prune(G: nx.Graph, seed: Optional[RandomState]) -> Generator[str, None, None]:
    node_to_remove = (
        pd.Series(dict(G.degree))
        .pipe(lambda x: x[x == x.max()])  # select the subjects with the most relatives
        .sample(n=1, random_state=seed)  # randomly select one of these subjects
        .index[0]
    )
    G.remove_node(node_to_remove)  # remove this subject

    non_connected_nodes = list(nx.isolates(G))  # Find subjects with no remaining relatives
    G.remove_nodes_from(non_connected_nodes)  # Remove these isolated subjects from the graph

    yield f"{node_to_remove} {node_to_remove}"  # add node_to_remove to pruning list


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({"out_file": Path(snakemake.output[0])})  # type: ignore # noqa
        defaults.update(snakemake.params)  # type: ignore # noqa
        main(**defaults)
    else:
        app()
