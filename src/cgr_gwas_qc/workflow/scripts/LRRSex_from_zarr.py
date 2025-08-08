from pathlib import Path
from typing import List

import numcodecs
import numpy as np
import sgkit as sg
import typer
import xarray as xr
from dask.distributed import Client

# import cgr_gwas_qc.parsers.virtual_zarr as virtual_zarr


app = typer.Typer(add_completion=False)


def virtual_zarr(
    filenames: list[str], sample_chunksize: int = 100, variant_chunksize: int = 100000
) -> xr.Dataset:
    """Read a Zarr dataset containing sample-level GWAS QC data.

    Args:
        filenames (str): Path(s) to the Zarr dataset.

    Returns:
        xr.Dataset: The dataset containing sample-level GWAS QC data.
    """

    if len(filenames) == 1:
        ds = sg.load_dataset(filenames[0], storage_options={"mode": "r"})
        return ds
    elif len(filenames) > 1:
        call_variables = [
            "call_BAF",
            "call_GQ",
            "call_IGC",
            "call_LRR",
            "call_NORMX",
            "call_NORMY",
            "call_R",
            "call_THETA",
            "call_X",
            "call_Y",
            "call_genotype",
            "call_genotype_mask",
            "call_genotype_phased",
        ]

        def _fix_filter_encoding(ds):
            for v in ds:
                # Workaround for https://github.com/pydata/xarray/issues/4380
                ds[v].encoding.pop("chunks", None)

                # Remove VLenUTF8 from filters to avoid double encoding error https://github.com/pydata/xarray/issues/3476
                filters = ds[v].encoding.get("filters", None)
                var_len_str_codec = numcodecs.VLenUTF8()
                if filters is not None and var_len_str_codec in filters:
                    filters = list(filters)
                    filters.remove(var_len_str_codec)
                    ds[v].encoding["filters"] = filters
            return ds

        filenames_meta = []
        for i in range(len(filenames)):
            ds_i = xr.open_dataset(filenames[i])
            filenames_meta.append(
                {
                    "path": filenames[i],
                    "n_variants": ds_i.sizes["variants"],
                    "n_samples": ds_i.sizes["samples"],
                }
            )
        ds = xr.open_mfdataset(
            filenames,
            concat_dim="samples",
            combine="nested",
            data_vars=call_variables + ["sample_id"],
            compat="override",
            coords="minimal",
            chunks={"samples": sample_chunksize},
        )
        ds = _fix_filter_encoding(ds)
        return ds


@app.command()
def main(
    zarr_ds: List[str],
    outfile: Path = typer.Option(
        ...,
        help="Path to output file to write sample median intensities.",
        file_okay=True,
        writable=True,
    ),
    threads: int = typer.Option(
        1,
        help="Number of CPU threads to use.",
    ),
    sd_cutoff: float = typer.Option(
        1.0,
        help="Standard deviation cutoff for determining sex.",
    ),
    median_cutoff: float = typer.Option(
        -1.6,
        help="Median LRR cutoff for determining sex.",
    ),
) -> None:

    Client(n_workers=threads, threads_per_worker=1)
    z = virtual_zarr(zarr_ds)

    Y_check = np.where(z.contig_id == "chrY")[0]

    if Y_check.size == 0:
        raise ValueError(
            "No 'chrY' contig found in the dataset. Please check the input Zarr dataset."
        )

    contig_id = Y_check[0]
    z = z.sel(variants=(z.variant_contig == contig_id).load())

    predicted_sex = np.empty(
        z.sample_id.shape,
        dtype=[
            ("sample_id", np.object_),
            ("chrY_sex", np.object_),
            ("median_LRR", np.float16),
            ("sd_LRR", np.float16),
        ],
    )
    predicted_sex["sample_id"] = z.sample_id.compute().values
    predicted_sex["median_LRR"] = z.call_LRR.median(axis=0).astype(np.float16).compute().values
    predicted_sex["sd_LRR"] = z.call_LRR.std(axis=0).astype(np.float16).compute().values
    predicted_sex["chrY_sex"] = xr.where(
        (predicted_sex["median_LRR"] < median_cutoff) | (predicted_sex["sd_LRR"] > sd_cutoff),
        "F",
        "M",
    )

    np.savetxt(
        outfile,
        predicted_sex,
        fmt=["%s", "%s", "%.3f", "%.3f"],
        delimiter=",",
        header="Sample_ID,chrY_sex,Median_LRR,SD_LRR",
        comments="",
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        main(
            **{k: v for k, v in snakemake.input.items() if not k.startswith("_")},  # type: ignore # noqa
            **{k: v for k, v in snakemake.params.items() if not k.startswith("_")},  # type: ignore # noqa
            outfile=snakemake.output[0],  # type: ignore # noqa
            threads=snakemake.threads,  # type: ignore # noqa
        )
    else:
        app()
