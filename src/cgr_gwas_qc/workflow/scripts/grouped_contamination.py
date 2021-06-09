import shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import List

import pandas as pd
from snakemake.rules import expand

from cgr_gwas_qc.typing import PathLike
from cgr_gwas_qc.workflow.scripts import agg_verifyidintensity, gtc2adpc, verifyidintensity


def main(
    sample_sheet_csv: PathLike,
    bpm_file: PathLike,
    abf_file: PathLike,
    grp: str,
    gtc_pattern: str,
    snps: int,
    conda_env: PathLike,
    outfile: PathLike,
    notemp: bool = False,
    threads: int = 8,
):
    ss = pd.read_csv(sample_sheet_csv).query(f"not is_missing_gtc and cluster_group == '{grp}'")

    # Make temp folders
    outdir = Path(outfile).parent

    tmp_adpc = outdir / "temp_adpc"
    tmp_adpc.mkdir(exist_ok=True, parents=True)

    tmp_contam = outdir / "temp_contam"
    tmp_contam.mkdir(exist_ok=True, parents=True)

    # Outfile patterns
    adpc_pattern = (tmp_adpc / "{Sample_ID}.adpc.bin").as_posix()
    contam_pattern = (tmp_contam / "{Sample_ID}.txt").as_posix()

    # Convert GTC to ADPC
    convert_gtc_to_adpc(ss, bpm_file, gtc_pattern, adpc_pattern, threads)

    # Estimate Contamination
    contam_files = estimate_contamination(
        ss, adpc_pattern, abf_file, contam_pattern, snps, conda_env, threads
    )
    agg_verifyidintensity.main(contam_files, outfile)

    if not notemp:
        shutil.rmtree(tmp_adpc)
        shutil.rmtree(tmp_contam)


def convert_gtc_to_adpc(ss, bpm_file, gtc_pattern, outfile_pattern, threads):
    with ProcessPoolExecutor(threads) as executor:
        futures = [
            executor.submit(
                gtc2adpc.main,
                gtc_pattern.format(**record.to_dict()),
                bpm_file,
                outfile_pattern.format(**record.to_dict()),
            )
            for _, record in ss.iterrows()
        ]

    # Wait for results
    [res.result() for res in futures]

    # Check the expected output files
    outfiles = [Path(f) for f in expand(outfile_pattern, **ss.to_dict("list"))]
    assert all(f.exists() for f in outfiles)


def estimate_contamination(
    ss, adpc_pattern, abf_file, outfile_pattern, snps, conda_env, threads
) -> List[Path]:
    with ProcessPoolExecutor(threads) as executor:
        futures = [
            executor.submit(
                verifyidintensity.main,
                adpc_pattern.format(**record.to_dict()),
                abf_file,
                outfile_pattern.format(**record.to_dict()),
                snps,
                conda_env,
            )
            for _, record in ss.iterrows()
        ]

    # Wait for results
    [res.result() for res in futures]

    # Check and Return the expected output files
    outfiles = [Path(f) for f in expand(outfile_pattern, **ss.to_dict("list"))]
    assert all(f.exists() for f in outfiles)
    return outfiles


if __name__ == "__main__" and "snakemake" in locals():
    main(
        **{k: v for k, v in snakemake.input.items() if not k.startswith("_")},  # type: ignore # noqa
        grp=snakemake.params.grp,  # type: ignore # noqa
        gtc_pattern=snakemake.params.gtc_pattern,  # type: ignore # noqa
        snps=int(snakemake.params.snps),  # type: ignore # noqa
        conda_env=snakemake.params.conda_env,  # type: ignore # noqa
        notemp=bool(snakemake.params.notemp),  # type: ignore # noqa
        outfile=snakemake.output[0],  # type: ignore # noqa
        threads=snakemake.threads or 2,  # type: ignore # noqa
    )
