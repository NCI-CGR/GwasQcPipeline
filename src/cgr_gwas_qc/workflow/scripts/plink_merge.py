import subprocess as sp
from pathlib import Path

import typer

app = typer.Typer(add_completion=False)


@app.command()
def main(merge_list: Path, out_prefix: str, threads: int = 8, mem_mb: int = 1024 * 8):
    """Merge PLINK PED/MAP files into a set of BED/BIM/FAM.

    Args:
        merge_list (Path): a list of files to merge in the format "{ped} {map}\n".
        out_prefix (str): the output name to use for BED/BIM/FAM/NOSEX/LOG files.
        threads (int, optional): Defaults to 8.
        mem_mb (int, optional): memory in MBs. Defaults to 1024*8.
    """
    cmd = [
        "plink",
        "--merge-list",
        merge_list.as_posix(),
        "--make-bed",
        "--out",
        out_prefix,
        "--threads",
        str(threads),
        "--memory",
        str(mem_mb),
    ]
    sp.run(cmd, check=True)


if __name__ == "__main__":
    if "snakemake" in locals():
        main(
            merge_list=Path(snakemake.input.merge_list),  # type: ignore # noqa
            out_prefix=snakemake.params.out_prefix,  # type: ignore # noqa
            threads=snakemake.threads,  # type: ignore # noqa
            mem_mb=snakemake.resources.mem_mb,  # type: ignore # noqa
        )
    else:
        app()
