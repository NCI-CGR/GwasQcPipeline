import subprocess as sp
from pathlib import Path
from typing import List

import typer

app = typer.Typer(add_completion=False)


@app.command()
def main(
    inputs: List[List[Path]],
    out_prefix: str,
    notemp: bool = False,
    threads: int = 8,
    mem_mb: int = 1024 * 8,
):
    """Merge PLINK PED/MAP files into a set of BED/BIM/FAM.

    Args:
        merge_list (Path): a list of files to merge in the format "{ped} {map}\n".
        out_prefix (str): the output name to use for BED/BIM/FAM/NOSEX/LOG files.
        threads (int, optional): Defaults to 8.
        mem_mb (int, optional): memory in MBs. Defaults to 1024*8.
    """
    outdir = Path(out_prefix).parent
    merge_list = create_merge_list(inputs, outdir)
    run_plink_merge(merge_list, out_prefix, threads, mem_mb)

    if not notemp:
        merge_list.unlink()


def create_merge_list(inputs: List[List[Path]], outdir: Path) -> Path:
    """Create a PLINK merge file.

    Creates a text file that plink uses when creating merged files. This script
    tries to capture all of PLINKs capabilities. PLINK will do slightly
    different things depending on how many file(s) are provided for each
    samples.

    - 1 file, is assumed to be the prefix for a binary fileset.
    - 2 files, are assumed to be the full filenames for a text fileset (.ped, .map).
    - 3 files, are assumed to be the full filenames for a binary fileset (.bed, .bim, .fam).

    Args:
        inputs (List[List[Path]]): A list of 1, 2 or 3 filesets.
        outfile (Path): File name to save the merge list.

    Reference:
        - https://www.cog-genomics.org/plink/1.9/data#merge_list
    """
    outfile = outdir / "plink_merge_list.txt"
    outfile.write_text(
        "\n".join(
            " ".join(filename.resolve().as_posix() for filename in fileset)
            for fileset in zip(*inputs)
        )
    )
    return outfile


def run_plink_merge(merge_list: Path, out_prefix: str, threads: int, mem_mb: int):
    """Run PLINK merge as a subprocess."""
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
    return sp.run(cmd, check=True)


if __name__ == "__main__":
    if "snakemake" in locals():
        if (input_items := list(snakemake.input.items())) :  # type: ignore # noqa
            # Here I am expecting a set of NamedLists.
            # The input directive must be in the format:
            # file1=[...],
            # file2=[...],
            # The names (file1, file2, ...) do not matter but PLINK does expect
            # a specific order, so refer to their documentation.

            # inputs = [[file1a, file1b, ...], [file2a, file2b,...], ...]
            inputs = [[Path(filename) for filename in named_list[1]] for named_list in input_items]

        else:
            # If just a list of files are provided, then plink assumes they are
            # file prefixes of BED/BIM/FAM

            # inputs = [[file1, file2, ...]]
            inputs = [[Path(filename) for filename in snakemake.input]]  # type: ignore # noqa

        main(
            inputs=inputs,
            out_prefix=snakemake.params.out_prefix,  # type: ignore # noqa
            notemp=snakemake.params.get("notemp", False),  # type: ignore # noqa
            threads=snakemake.threads,  # type: ignore # noqa
            mem_mb=snakemake.resources.mem_mb,  # type: ignore # noqa
        )
    else:
        app()
