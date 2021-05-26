import shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import List

import pandas as pd
from snakemake.rules import expand

from cgr_gwas_qc.typing import PathLike
from cgr_gwas_qc.workflow.scripts import agg_median_idat_intensity, median_idat_intensity


def main(
    sample_sheet_csv: PathLike,
    grp: str,
    idat_red_pattern: str,
    idat_green_pattern: str,
    rscript: PathLike,
    conda_env: PathLike,
    outfile: PathLike,
    notemp: bool = False,
    threads: int = 8,
):
    ss = pd.read_csv(sample_sheet_csv).query(f"cluster_group == '{grp}'")
    tmp_path = setup_folders(outfile)

    # Calculate median IDAT intensity
    intensity_pattern = (tmp_path / "{Sample_ID}.csv").as_posix()
    intensity_files = calculate_median_idat_intensity(
        ss, idat_red_pattern, idat_green_pattern, intensity_pattern, rscript, conda_env, threads
    )
    agg_median_idat_intensity.main(intensity_files, outfile)

    if not notemp:
        shutil.rmtree(tmp_path)


def setup_folders(outfile: PathLike) -> Path:
    """Set-up my output folder"""
    outdir = Path(outfile).parent
    tmp_path = outdir / "temp/median_idat"
    tmp_path.mkdir(exist_ok=True, parents=True)
    return tmp_path


def calculate_median_idat_intensity(
    ss, red_pattern, green_pattern, outfile_pattern, rscript, conda_env, threads
) -> List[Path]:
    with ProcessPoolExecutor(threads) as executor:
        futures = [
            executor.submit(
                median_idat_intensity.main,
                red_pattern.format(**record.to_dict()),
                green_pattern.format(**record.to_dict()),
                record.Sample_ID,
                rscript,
                conda_env,
                outfile_pattern.format(**record.to_dict()),
            )
            for _, record in ss.iterrows()
        ]

    # Wait for results
    [res.result() for res in futures]

    # Check and Return files
    outfiles = [Path(f) for f in expand(outfile_pattern, **ss.to_dict("list"))]
    assert all(f.exists() for f in outfiles)
    return outfiles


if __name__ == "__main__" and "snakemake" in locals():
    main(
        sample_sheet_csv=snakemake.input.sample_sheet_csv,  # type: ignore # noqa
        **{k: v for k, v in snakemake.params.items()},  # type: ignore # noqa
        outfile=snakemake.output[0],  # type: ignore # noqa
        threads=snakemake.threads or 2,  # type: ignore # noqa
    )
