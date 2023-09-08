import subprocess as sp
from pathlib import Path

from cgr_gwas_qc.conda import get_snakemake_conda_env
from cgr_gwas_qc.typing import PathLike


def main(
    red: PathLike,
    green: PathLike,
    sample_id: str,
    rscript: PathLike,
    conda_env: PathLike,
    outfile: PathLike,
):
    Rscript = Path(get_snakemake_conda_env(conda_env)) / "bin/Rscript"
    cmd = f"{Rscript} --vanilla {rscript} {sample_id} {red} {green} {outfile}"
    sp.check_output(cmd, shell=True, executable="/bin/bash")


if __name__ == "__main__" and "snakemake" in locals():
    main(
        **{k: v for k, v in snakemake.input.items() if not k.startswith("_")},  # type: ignore # noqa
        **{k: v for k, v in snakemake.params.items()},  # type: ignore # noqa
        outfile=snakemake.output[0],  # type: ignore # noqa
    )
