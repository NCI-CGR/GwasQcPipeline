import subprocess as sp

from cgr_gwas_qc.conda import conda_activate
from cgr_gwas_qc.typing import PathLike


def main(
    red: PathLike,
    green: PathLike,
    sample_id: str,
    rscript: PathLike,
    conda_env: PathLike,
    outfile: PathLike,
):
    cmd = (
        f"{conda_activate(conda_env)}"
        "&& Rscript --vanilla"
        f" {rscript}"
        f" {sample_id}"
        f" {red}"
        f" {green}"
        f" {outfile}"
    )
    sp.check_output(cmd, shell=True)


if __name__ == "__main__" and "snakemake" in locals():
    main(
        **{k: v for k, v in snakemake.input.items() if not k.startswith("_")},  # type: ignore # noqa
        **{k: v for k, v in snakemake.params.items()},  # type: ignore # noqa
        outfile=snakemake.output[0],  # type: ignore # noqa
    )
