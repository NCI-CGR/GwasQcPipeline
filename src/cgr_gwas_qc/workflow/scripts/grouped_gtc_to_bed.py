import shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import List

import pandas as pd
from snakemake.rules import expand

from cgr_gwas_qc.typing import PathLike
from cgr_gwas_qc.workflow.scripts import gtc2plink, plink_merge


def main(
    sample_sheet_csv: PathLike,
    bpm_file: PathLike,
    grp: str,
    gtc_pattern: str,
    strand: str,
    out_prefix: str,
    conda_env: PathLike,
    notemp: bool = False,
    threads: int = 8,
    mem_mb: int = 1024 * 8,
):
    """Grouped version of GTC to BED/BIM/FAM conversion."""
    ss = pd.read_csv(sample_sheet_csv).query(f"cluster_group == '{grp}'")
    tmp_dir = Path(out_prefix).parent / "temp_gtc2bed"
    tmp_dir.mkdir(exist_ok=True, parents=True)

    ped_map_files = convert_gtc_to_ped_map(ss, bpm_file, strand, tmp_dir, gtc_pattern, threads)
    plink_merge.main(ped_map_files, out_prefix, conda_env, notemp, threads, mem_mb)

    # Clean-up
    if not notemp:
        shutil.rmtree(tmp_dir)  # remove per sample ped/map files and merge list
        Path(out_prefix).with_suffix(".log").unlink()  # remove the merge log


def convert_gtc_to_ped_map(ss, bpm_file, strand, outdir, gtc_pattern, threads) -> List[List[Path]]:
    with ProcessPoolExecutor(threads) as executor:
        futures = [
            executor.submit(
                gtc2plink.main,
                gtc_pattern.format(**record.to_dict()),
                bpm_file,
                record.Sample_ID,
                outdir,
                strand,
                None,
            )
            for _, record in ss.iterrows()
        ]

    # Wait for results
    [res.result() for res in futures]

    # Check and return results
    ped = [Path(f) for f in expand(outdir / "{Sample_ID}.ped", **ss.to_dict("list"))]
    map_ = [Path(f) for f in expand(outdir / "{Sample_ID}.map", **ss.to_dict("list"))]

    assert all(f.exists() for f in ped)
    assert all(f.exists() for f in map_)

    return [ped, map_]


if __name__ == "__main__" and "snakemake" in locals():
    main(
        **{k: v for k, v in snakemake.input.items() if not k.startswith("_")},  # type: ignore # noqa
        **{k: v for k, v in snakemake.params.items()},  # type: ignore # noqa
        threads=snakemake.threads or 2,  # type: ignore # noqa
        mem_mb=snakemake.resources.get("mem_mb", 1024 * 2),  # type: ignore # noqa
    )
