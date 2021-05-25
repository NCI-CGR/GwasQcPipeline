import shutil
import subprocess as sp
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import List

import pandas as pd
import typer

from cgr_gwas_qc.conda import conda_activate
from cgr_gwas_qc.workflow.scripts import agg_median_idat_intensity

app = typer.Typer(add_completion=False)

SCRIPT_DIR = Path(__file__).resolve().parent


@app.command()
def main(
    sample_sheet_csv: Path,
    grp: str,
    idat_red_pattern: str,
    idat_green_pattern: str,
    conda_env: str,
    outfile: Path,
    notemp: bool = False,
    threads: int = 8,
):
    tmp_path = setup_folders(outfile)
    intensity_files = calculate_median_idat_intensity(
        sample_sheet_csv, grp, idat_red_pattern, idat_green_pattern, conda_env, tmp_path, threads
    )
    agg_median_idat_intensity.main(intensity_files, outfile)

    if not notemp:
        shutil.rmtree(tmp_path)


def setup_folders(outfile: Path) -> Path:
    """Set-up my output folder"""
    outdir = outfile.parent
    tmp_path = outdir / "temp"
    tmp_path.mkdir(exist_ok=True, parents=True)
    return tmp_path


def calculate_median_idat_intensity(
    sample_sheet_csv, grp, red_pattern, green_pattern, conda_env, outdir, threads
) -> List[Path]:
    ss = pd.read_csv(sample_sheet_csv, index_col=["Sample_ID"]).query(f"cluster_group == '{grp}'")
    with ProcessPoolExecutor(threads) as executor:
        futures = [
            executor.submit(
                _run_illuminaio,
                record,
                red_pattern,
                green_pattern,
                conda_env,
                outdir,
            )
            for _, record in ss.iterrows()
        ]
    return [res.result() for res in futures]


def _run_illuminaio(record, red_pattern, green_pattern, conda_env, outdir) -> Path:
    """Run the R package illuminaio

    The run command looks a little confusing b/c I have to activate the correct
    conda environment first.
    """
    sample_id = record.name

    if "SentrixBarcode_A" in record and "SentrixPosition_A" in record:
        chip_id = "_".join([str(record.SentrixBarcode_A), str(record.SentrixPosition_A)])
    else:
        chip_id = sample_id

    red_file = Path(red_pattern.format(**record.to_dict()))
    green_file = Path(green_pattern.format(**record.to_dict()))
    outfile = outdir / f"{sample_id}.csv"

    cmd = (
        f"{conda_activate(conda_env)}"
        "&& Rscript --vanilla"
        f" {SCRIPT_DIR / 'median_idat_intensity.R'}"
        f" {sample_id}"
        f" {chip_id}"
        f" {red_file}"
        f" {green_file}"
        f" {outfile}"
    )
    sp.check_output(cmd, shell=True)

    return outfile


if __name__ == "__main__":
    if "snakemake" in locals():
        main(
            sample_sheet_csv=Path(snakemake.input.sample_sheet_csv),  # type: ignore # noqa
            **{k: v for k, v in snakemake.params.items()},  # type: ignore # noqa
            outfile=Path(snakemake.output[0]),  # type: ignore # noqa
            threads=snakemake.threads or 2,  # type: ignore # noqa
        )
    else:
        app()
