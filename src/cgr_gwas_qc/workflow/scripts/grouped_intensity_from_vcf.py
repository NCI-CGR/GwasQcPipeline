#!/usr/bin/env python
# coding: utf-8


import shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import List

from snakemake.rules import expand

from cgr_gwas_qc.parsers import sample_sheet
from cgr_gwas_qc.typing import PathLike
from cgr_gwas_qc.workflow.scripts import agg_median_idat_intensity, median_intensity_from_vcf


def main(
    vcf_file: PathLike,
    sample_sheet_csv: PathLike,
    grp: str,
    outfile: PathLike,
    notemp: bool = False,
    threads: int = 8,
):
    ss = sample_sheet.read(sample_sheet_csv).query(f"cluster_group == '{grp}'")
    tmp_dir = Path(outfile).parent / "temp_median_idat"
    tmp_dir.mkdir(exist_ok=True, parents=True)

    intensity_files = calculate_median_intensity_from_vcf(ss, vcf_file, tmp_dir, threads)

    agg_median_idat_intensity.main(intensity_files, outfile)

    if not notemp:
        shutil.rmtree(tmp_dir)


def calculate_median_intensity_from_vcf(ss, vcf_file, outdir, threads) -> List[Path]:
    outfile_pattern = (outdir / "{Sample_ID}.csv").as_posix()
    with ProcessPoolExecutor(threads) as executor:
        futures = [
            executor.submit(
                median_intensity_from_vcf.main,
                vcf_file,
                record.Sample_ID,
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
        **{k: v for k, v in snakemake.input.items()},  # type: ignore # noqa
        **{k: v for k, v in snakemake.params.items()},  # type: ignore # noqa
        outfile=snakemake.output[0],  # type: ignore # noqa
        threads=snakemake.threads or 2,  # type: ignore # noqa
    )
