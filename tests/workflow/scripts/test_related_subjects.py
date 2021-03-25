import networkx as nx
import numpy as np

from cgr_gwas_qc.workflow.scripts.related_subjects import create_prune_list, create_qc_families


def test_create_prune_list():
    pairwise_related_subjects = [
        ("S001", "S002"),
        ("S001", "S003"),
        ("S001", "S004"),
        ("S005", "S006"),
        ("S005", "S007"),
        ("S008", "S009"),
        ("S010", "S011"),
    ]

    G = nx.Graph()
    G.add_edges_from(pairwise_related_subjects)

    prune_list = list(create_prune_list(G, seed=np.random.RandomState(42)))
    assert sorted(["S001 S001", "S005 S005", "S009 S009", "S011 S011"]) == sorted(prune_list)


def test_summarize_relatives():
    # GIVEN: a set of disconnected relatives
    pairwise_related_subjects = [
        ("S001", "S002"),
        ("S001", "S003"),
        ("S001", "S004"),
        ("S005", "S006"),
        ("S005", "S007"),
        ("S008", "S009"),
        ("S010", "S011"),
    ]

    G = nx.Graph()
    G.add_edges_from(pairwise_related_subjects)

    # WHEN:
    obs = create_qc_families(G)
    exp_ = ["S001|S002|S003|S004", "S005|S006|S007", "S008|S009", "S010|S011"]
    assert sorted(exp_) == sorted(obs.values)
