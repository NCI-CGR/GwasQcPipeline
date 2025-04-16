#!/usr/bin/env python3
import concurrent.futures
import subprocess as sp
import tempfile
from functools import partial
from io import StringIO
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import sgkit as sg
import typer
from dask.distributed import Client
from numba import njit

app = typer.Typer(add_completion=True)


@app.command()
def main(
    zarr_ds: Path = typer.Option(
        ...,
        help="Path to a multisample zarr dataset converted from VCF/BCF file.",
        exists=True,
        readable=True,
    ),
    outfile: Path = typer.Option(
        ...,
        help="Path to output file to write ADPC.bin intensities.",
        file_okay=True,
        writable=True,
    ),
    abf: Optional[Path] = typer.Option(
        None,
        help="Path to abf.txt file containing B-allele frequencies. If None, conamianation check will be performed without it.",
        exists=True,
        readable=True,
    ),
    adpc: Optional[Path] = typer.Option(
        None,
        help="Path to output ADPC.bin file. If None, a temporary file will be created.",
        exists=True,
        readable=True,
    ),
    batch_size: int = typer.Option(
        100,
        help="Number of samples to process in a batch. More samples per batch will use more memory.",
    ),
    threads: int = typer.Argument(
        1,
        help="Number of threads to use for parallel processing.",
    ),
) -> None:

    def check_verifyidintensity(adpc: Path, n_markers: int, sampleidx: int, abf: Path = None):
        """Checks contamination for a single sample from a multisample ADPC.bin file using verifyIDintensity tool.
        Parameters
        ----------
        adpc : Path
            Path to ADPC.bin file.
        n_markers : int
            Number of markers in the ADPC.bin file.
        sampleidx : int
            Index of the sample to check for contamination.
        abf : Path
            Path to abf.txt file containing B-allele frequencies. If None, conamianation check will be performed without it.
        Returns
        -------
        pd.DataFrame
            Contamination check results.
        """
        s_blocksize = 18 * n_markers
        skip = 16 + s_blocksize * sampleidx
        with tempfile.NamedTemporaryFile(delete=True) as temp_sample_adpc:
            cmd_add_header = f"dd skip=0 count=16 iflag=skip_bytes,count_bytes if={adpc} > {temp_sample_adpc.name}"
            cmd_retrieve_sample = f"dd skip={skip} count={s_blocksize} iflag=skip_bytes,count_bytes if={adpc} >> {temp_sample_adpc.name}"
            if abf is not None:
                cmd_verifyidintensity_check = f"verifyIDintensity -n 1 -v -p -m {n_markers} -i {temp_sample_adpc.name} -b {abf}"
                rows_to_skip = [0, 1, 3]
            else:
                cmd_verifyidintensity_check = (
                    f"verifyIDintensity -n 1 -v -p -m {n_markers} -i {temp_sample_adpc.name}"
                )
                rows_to_skip = [1]
            cmd = cmd_add_header + ";" + cmd_retrieve_sample + ";" + cmd_verifyidintensity_check
            contam_out = sp.run(cmd, shell=True, capture_output=True, text=True, check=True)
            return pd.read_csv(StringIO(contam_out.stdout), sep="\\s+", skiprows=rows_to_skip)

    def zarr2adpc(zarr_ds: Path, adpc: Path, batch_size: int, threads: int):
        """Converts a multisample zarr dataset to ADPC.bin file.
        Parameters
        ----------
        zarr_ds : Path
            Path to a multisample zarr dataset converted from VCF/BCF file.
        adpc : Path
            Path to output ADPC.bin file.
        batch_size : int
            Number of samples to process in a batch.
        threads : int
            Number of threads to use for parallel processing.
        Returns
        -------
        Path
            Path to ADPC.bin file.
        np.ndarray
            Sample IDs.
        int
            Number of markers in the ADPC.bin file.
        """

        def get_illumina_genotype(allele_a, allele_b, vcf_gt):
            illumina_genotype = np.full(vcf_gt.shape[0], 3, dtype=np.uint16)
            for v in np.ndindex(illumina_genotype.shape[0]):
                if vcf_gt[v][0] != vcf_gt[v][1]:
                    illumina_genotype[v] = 1
                elif vcf_gt[v][0] == vcf_gt[v][1] == allele_a[v]:
                    illumina_genotype[v] = 0
                elif vcf_gt[v][0] == vcf_gt[v][1] == allele_b[v]:
                    illumina_genotype[v] = 2
            return illumina_genotype

        def get_illumina_genotype_multisample(allele_a, allele_b, vcf_gt):
            illumina_genotype = np.full(vcf_gt.shape[:2], 3, dtype=np.uint16)
            for s in np.ndindex(illumina_genotype.shape[1]):
                for v in np.ndindex(illumina_genotype.shape[0]):
                    if vcf_gt[v][s][0] != vcf_gt[v][s][1]:
                        illumina_genotype[v][s] = 1
                    elif vcf_gt[v][s][0] == vcf_gt[v][s][1] == allele_a[v]:
                        illumina_genotype[v][s] = 0
                    elif vcf_gt[v][s][0] == vcf_gt[v][s][1] == allele_b[v]:
                        illumina_genotype[v][s] = 2
            return illumina_genotype

        record_struct = {
            "X": "<H",
            "Y": "<H",
            "NORMX": "<f",
            "NORMY": "<f",
            "IGC": "<f",
            "illumina_GT": "<H",
        }
        Path(adpc).parent.mkdir(parents=True, exist_ok=True)
        with open(adpc, "wb") as f:
            f.write(np.repeat(0, 8).astype("<H").tobytes())
        Client(n_workers=threads, threads_per_worker=1)
        variables = [
            "call_X",
            "call_Y",
            "call_NORMX",
            "call_NORMY",
            "call_IGC",
            "variant_ALLELE_A",
            "variant_ALLELE_B",
            "call_genotype",
        ]
        z = sg.load_dataset(zarr_ds, storage_options={"mode": "r"})
        sample_ids = z.sample_id.to_numpy()
        n_markers = z.variants.shape[0]
        z["call_X"] = z.call_X.astype(np.uint16)
        z["call_Y"] = z.call_Y.astype(np.uint16)
        z["variant_ALLELE_A"] = z.variant_ALLELE_A.astype(np.int8)
        z["variant_ALLELE_B"] = z.variant_ALLELE_B.astype(np.int8)
        z = z[variables]
        with open(adpc, "ab") as f:
            for i in range(0, z.samples.shape[0], batch_size):
                z_batch = z.isel(samples=slice(i, i + batch_size))
                z_batch["call_illumina_genotype"] = (
                    ("variants", "samples"),
                    njit(get_illumina_genotype_multisample, parallel=False)(
                        z_batch.variant_ALLELE_A.values,
                        z_batch.variant_ALLELE_B.values,
                        z_batch.call_genotype.values,
                    ),
                )
                z_batch = z_batch.drop_vars(
                    ["call_genotype", "variant_ALLELE_A", "variant_ALLELE_B"]
                )
                z_batch = z_batch.load()
                z_batch["call_NORMX"] = (
                    ("variants", "samples"),
                    np.nan_to_num(z_batch.call_NORMX, posinf=0, neginf=0, copy=False),
                )
                z_batch["call_NORMY"] = (
                    ("variants", "samples"),
                    np.nan_to_num(z_batch.call_NORMY, posinf=0, neginf=0, copy=False),
                )
                f.write(
                    z_batch.to_dataframe()
                    .sort_index(level=["samples", "variants"])
                    .to_records(index=False, column_dtypes=record_struct)
                    .tobytes()
                )
                f.flush()
        return adpc, sample_ids, n_markers

    if adpc is None:
        adpc_temp = tempfile.NamedTemporaryFile(delete=False)
        adpc = Path(adpc_temp.name)
    adpc, sample_ids, n_markers = zarr2adpc(zarr_ds, adpc, batch_size, threads)
    check_verifyidintensity = partial(check_verifyidintensity, adpc, n_markers, abf=abf)
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        results = executor.map(check_verifyidintensity, range(0, sample_ids.shape[0]))
    results = pd.concat(results).assign(ID=sample_ids).rename({"ID": "Sample_ID"}, axis=1)
    results.to_csv(outfile, index=False)
    if adpc_temp is None:
        adpc_temp.unlink()


if __name__ == "__main__":
    if "snakemake" in locals():
        defaults = {}
        defaults.update({k: Path(v) for k, v in snakemake.input.items()})  # type: ignore # noqa
        defaults.update({k: v for k, v in snakemake.params.items()})  # type: ignore # noqa
        defaults.update({"outfile": Path(snakemake.output[0])})  # type: ignore # noqa
        defaults.update({"threads": snakemake.threads})  # type: ignore # noqa
        main(**defaults)
    else:
        app()
