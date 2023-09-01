import subprocess as sp
from pathlib import Path
from typing import Sequence

from cgr_gwas_qc.conda import conda_activate
from cgr_gwas_qc.typing import PathLike

NestedPaths = Sequence[Sequence[PathLike]]


def main(
    inputs: NestedPaths,
    out_prefix: str,
    conda_env: PathLike,
    notemp: bool = False,
    threads: int = 8,
    mem_mb: int = 1024 * 8,
):
    """Merge PLINK PED/MAP files into a set of BED/BIM/FAM.

    Args:
        merge_list (Path): a list of files to merge in the format "{ped} {map}\n".
        out_prefix (str): the output name to use for BED/BIM/FAM/NOSEX/LOG files.
        conda_env (PathLike): The path to the conda env with PLINK installed.
        threads (int, optional): Defaults to 8.
        mem_mb (int, optional): memory in MBs. Defaults to 1024*8.
    """
    outdir = Path(out_prefix).parent
    merge_list = create_merge_list(inputs, outdir)
    run_plink_merge(merge_list, out_prefix, conda_env, threads, mem_mb)

    if not notemp:
        merge_list.unlink()


def create_merge_list(inputs: NestedPaths, outdir: Path) -> Path:
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
            " ".join(Path(filename).resolve().as_posix() for filename in fileset)
            for fileset in zip(*inputs)
        )
    )
    return outfile


def run_plink_merge(merge_list, out_prefix, conda_env, threads, mem_mb):
    """Run PLINK merge as a subprocess."""
    cmd = (
        f"{conda_activate(conda_env)}"
        " && plink"
        " --make-bed"
        f" --merge-list {merge_list}"
        f" --out {out_prefix}"
        f" --threads {threads}"
        f" --memory {mem_mb}"
    )
    return sp.check_output(cmd, shell=True, executable="/bin/bash")


if __name__ == "__main__" and "snakemake" in locals():
    main(
        inputs=[v for k, v in snakemake.input.items() if not k.startswith("_")],  # type: ignore # noqa
        out_prefix=snakemake.params.out_prefix,  # type: ignore # noqa
        conda_env=snakemake.params.conda_env,  # type: ignore # noqa
        notemp=snakemake.params.get("notemp", False),  # type: ignore # noqa
        threads=snakemake.threads,  # type: ignore # noqa
        mem_mb=snakemake.resources.mem_mb,  # type: ignore # noqa
    )
