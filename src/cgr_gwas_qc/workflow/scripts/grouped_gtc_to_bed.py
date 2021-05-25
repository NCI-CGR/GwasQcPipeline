import shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import List, Tuple

import pandas as pd
import typer
from more_itertools import unzip

from cgr_gwas_qc.workflow.scripts import gtc2plink, plink_merge

app = typer.Typer(add_completion=False)


@app.command()
def main(
    sample_sheet_csv: Path,
    bpm_file: Path,
    grp: str,
    gtc_pattern: str,
    strand: str,
    out_prefix: str,
    notemp: bool = False,
    threads: int = 8,
    mem_mb: int = 1024 * 8,
):
    """Grouped version of GTC to BED/BIM/FAM conversion."""
    tmp_path = setup_folders(out_prefix)

    ped_map_files = convert_gtc_to_ped_map(
        sample_sheet_csv, bpm_file, grp, strand, tmp_path, gtc_pattern, threads
    )
    plink_merge.main(ped_map_files, out_prefix, threads, mem_mb)

    # Clean-up
    if not notemp:
        shutil.rmtree(tmp_path)  # remove per sample ped/map files and merge list


def setup_folders(out_prefix: str) -> Path:
    """Set-up my output folder"""
    out_dir = Path(out_prefix).parent
    tmp_path = out_dir / "temp"
    tmp_path.mkdir(exist_ok=True, parents=True)
    return tmp_path


def _gtc_to_ped_map(record, bpm_file, strand, out_dir, gtc_pattern) -> Tuple[Path, Path]:
    sample_id = record.name
    gtc_file = Path(gtc_pattern.format(**record.to_dict()))
    gtc2plink.main(gtc_file, bpm_file, sample_id, out_dir, strand, None)
    return out_dir / f"{sample_id}.ped", out_dir / f"{sample_id}.map"


def convert_gtc_to_ped_map(
    sample_sheet: Path,
    bpm_file: Path,
    grp: str,
    strand: str,
    out_dir: Path,
    gtc_pattern: str,
    threads: int,
) -> List[List[Path]]:
    ss = pd.read_csv(sample_sheet, index_col=["Sample_ID"]).query(f"cluster_group == '{grp}'")
    with ProcessPoolExecutor(threads) as executor:
        futures = [
            executor.submit(
                _gtc_to_ped_map,
                record,
                bpm_file,
                strand,
                out_dir,
                gtc_pattern,
            )
            for _, record in ss.iterrows()
        ]
    return list(map(list, unzip([res.result() for res in futures])))


if __name__ == "__main__":
    if "snakemake" in locals():
        main(
            sample_sheet_csv=Path(snakemake.input.sample_sheet_csv),  # type: ignore # noqa
            bpm_file=Path(snakemake.input.bpm_file),  # type: ignore # noqa
            grp=snakemake.params.grp,  # type: ignore # noqa
            gtc_pattern=snakemake.params.gtc_pattern,  # type: ignore # noqa
            strand=snakemake.params.strand,  # type: ignore # noqa
            out_prefix=snakemake.params.out_prefix,  # type: ignore # noqa
            notemp=snakemake.params.get("notemp", False),  # type: ignore # noqa
            threads=snakemake.threads or 2,  # type: ignore # noqa
            mem_mb=snakemake.resources.get("mem_mb", 1024 * 2),  # type: ignore # noqa
        )
    else:
        app()
