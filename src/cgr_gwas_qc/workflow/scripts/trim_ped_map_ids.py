import re
from typing import Generator

from cgr_gwas_qc.typing import PathLike


def main(ped: PathLike, map_: PathLike, out_ped: PathLike, out_map: PathLike, id_length: int = 39):
    with open(out_ped, "w") as fh:
        for row in trim_ped_ids(ped, id_length):
            fh.write(row)

    with open(out_map, "w") as fh:
        for row in trim_map_ids(map_, id_length):
            fh.write(row)


def trim_ped_ids(filename: PathLike, id_length: int) -> Generator[str, None, None]:
    with open(filename) as fh:
        for row in fh:
            yield _trim_ped_row(row, id_length)


def trim_map_ids(filename: PathLike, id_length: int) -> Generator[str, None, None]:
    with open(filename) as fh:
        for row in fh:
            yield _trim_map_row(row, id_length)


def _trim_id(id_: str, id_length: int) -> str:
    if len(id_) > id_length:
        return id_[:id_length]
    return id_


def _trim_ped_row(row: str, id_length: int) -> str:
    columns = re.split(r"[\t\s+]", row.strip())
    columns[0] = _trim_id(columns[0], id_length)
    columns[1] = _trim_id(columns[1], id_length)
    return " ".join(columns) + "\n"


def _trim_map_row(row: str, id_length: int) -> str:
    columns = re.split(r"[\t\s+]", row.strip())
    columns[1] = _trim_id(columns[1], id_length)
    return "\t".join(columns) + "\n"


if __name__ == "__main__" and "snakemake" in locals():
    main(
        snakemake.input.ped,  # type: ignore # noqa
        snakemake.input.map_,  # type: ignore # noqa
        snakemake.output.ped,  # type: ignore # noqa
        snakemake.output.map_,  # type: ignore # noqa
    )
