import subprocess as sp

from cgr_gwas_qc.conda import conda_activate
from cgr_gwas_qc.typing import PathLike


def main(adpc: PathLike, abf: PathLike, outfile: PathLike, snps: int, conda_env: PathLike):
    cmd = (
        f"{conda_activate(conda_env)}"
        " && verifyIDintensity -n 1 -v -p"
        f" -m {snps}"
        f" -b {abf}"
        f" -i {adpc}"
        f" > {outfile}"
    )
    sp.check_output(cmd, shell=True)


if __name__ == "__main__" and "snakemake" in locals():
    main(
        **{k: v for k, v in snakemake.input.items() if not k.startswith("_")},  # type: ignore # noqa
        **{k: int(v) if k == "snps" else v for k, v in snakemake.params.items()},  # type: ignore # noqa
        outfile=snakemake.output[0],  # type: ignore # noqa
    )
