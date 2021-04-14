import pandas as pd

from cgr_gwas_qc.typing import PathLike


def infer_relationship(king_coefficient: float):
    if king_coefficient > 0.354:
        return "ID"  # identical or MZ twin

    if 0.177 < king_coefficient <= 0.354:
        return "D1"  # 1st degree relative

    if 0.0884 < king_coefficient <= 0.177:
        return "D2"  # 2nd degree relative

    if 0.0442 < king_coefficient <= 0.0884:
        return "D3"  # 3rd degree relative

    return "UN"  # unrelated


def read_kinship(filename: PathLike) -> pd.DataFrame:
    dtypes = {
        "ID1": "string",
        "ID2": "string",
        "N_SNP": "UInt32",
        "HetHet": "float",
        "IBS0": "float",
        "Kinship": "float",
        "relationship": "string",
    }
    return (
        pd.read_csv(filename, sep="\t")
        .assign(relationship=lambda x: x.Kinship.apply(infer_relationship))
        .reindex(dtypes.keys(), axis=1)
    )
