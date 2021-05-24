from pathlib import Path
from typing import List

from typer import Typer

app = Typer(add_completion=False)


@app.command()
def main(inputs: List[List[Path]], outfile: Path):
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
    outfile.write_text(
        "\n".join(
            " ".join(filename.resolve().as_posix() for filename in fileset)
            for fileset in zip(*inputs)
        )
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        if (input_items := list(snakemake.input.items())) :  # type: ignore # noqa
            # Here I am expecting a set of NamedLists.
            # The inpuyt directive you must be in the format:
            # file1=[...],
            # file2=[...],
            # The names (file1, file2, ...) do not matter but PLINK does expect
            # a specific order, so refer to their documentation.

            # inputs = [[file1a, file1b, ...], [file2a, file2b,...], ...]
            inputs = [[Path(filename) for filename in named_list[1]] for named_list in input_items]

        else:
            # Here I am expected just a list of file names which plink will
            # assume are file prefixes of BED/BIM/FAM

            # inputs = [[file1, file2, ...]]
            inputs = [[Path(filename) for filename in snakemake.input]]  # type: ignore # noqa

        main(
            inputs=inputs,
            outfile=Path(snakemake.output[0]),  # type: ignore # noqa
        )
    else:
        app()
