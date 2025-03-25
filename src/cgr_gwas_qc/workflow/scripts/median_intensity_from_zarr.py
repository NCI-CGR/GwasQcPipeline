from pathlib import Path

import numpy as np
import sgkit as sg
import typer
import xarray as xr
from dask.distributed import Client

app = typer.Typer(add_completion=False)


@app.command()
def main(
    zarr_ds: Path = typer.Argument(
        ...,
        help="Path to a multisample zarr dataset converted from VCF/BCF file.",
        exists=True,
        readable=True,
    ),
    outfile: Path = typer.Argument(
        ...,
        help="Path to output file to write sample median intensities.",
        file_okay=True,
        writable=True,
    ),
    threads: int = typer.Argument(
        1,
        help="Number of CPU threads to use.",
    ),
) -> None:

    Client(n_workers=threads, threads_per_worker=1)
    z = sg.load_dataset(zarr_ds, storage_options={"mode": "r"})
    z["call_X"] = z.call_X.astype(np.uint16)
    z["call_Y"] = z.call_Y.astype(np.uint16)
    z = z.assign(xy_sum=z.call_X + z.call_Y)
    median_intensities = np.empty(
        z.sample_id.shape, dtype=[("Sample_ID", np.object_), ("median_intensity", np.float16)]
    )
    median_intensities["Sample_ID"] = z.sample_id
    median_intensities["median_intensity"] = (
        xr.DataArray.median(z.call_X + z.call_Y, dim="variants").compute().values
    )
    np.savetxt(
        outfile,
        median_intensities,
        fmt=["%s", "%.1f"],
        delimiter=",",
        header="Sample_ID,median_intensity",
        comments="",
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        main(
            **{k: v for k, v in snakemake.input.items() if not k.startswith("_")},  # type: ignore # noqa
            outfile=snakemake.output[0],  # type: ignore # noqa
            threads=snakemake.threads,  # type: ignore # noqa
        )
    else:
        app()
