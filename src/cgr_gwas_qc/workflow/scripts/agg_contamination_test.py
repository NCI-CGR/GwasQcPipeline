import re
from io import StringIO
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import typer

app = typer.Typer(add_completion=False)


@app.command()
def main(
    contamination_files: List[Path] = typer.Argument(
        ..., help="List of Paths with sample contamination estimates."
    ),
    median_intensity_file: Path = typer.Argument(
        ..., help="Path to aggregated median IDAT intensities.", exists=True, readable=True
    ),
    imiss_file: Path = typer.Argument(
        ...,
        help="Path to sample based missing report from lvl2 call rate.",
        exists=True,
        readable=True,
    ),
    intensity_threshold: int = typer.Argument(..., help="IDAT intensity threshold."),
    out_file: Path = typer.Argument(
        ..., help="Path to output file.", file_okay=True, writable=True
    ),
) -> None:
    contamination = pd.concat(
        [parse_sample_contamination_output(sample) for sample in contamination_files],
        ignore_index=True,
    )
    median_intensity = parse_median_intensity(median_intensity_file)
    imiss = parse_imiss(imiss_file)

    df = contamination.merge(median_intensity, on="Sample_ID").merge(imiss, on="Sample_ID")

    # Set %Mix to NA if < threshold or not in imiss
    mask = (df.median_intensity < intensity_threshold) & df.F_MISS.isnull()
    df.loc[mask, "%Mix"] = np.nan

    df[["Sample_ID", "%Mix", "LLK", "LLK0"]].to_csv(out_file, index=False)


def parse_sample_contamination_output(file_name: Path) -> pd.DataFrame:
    sample_id = file_name.stem.split(".")[0]
    data = remove_header_delim(file_name.read_text())
    contam = pd.read_csv(StringIO(data), skiprows=2, delim_whitespace=True)
    contam.ID = sample_id
    return contam.rename({"ID": "Sample_ID"}, axis=1)


def remove_header_delim(data: str) -> str:
    """Remove the table header delimiter from verifyIDintensity output.

    verifyIDintensity adds a row of "---" to deliminate the header. This
    makes it harder to read it into a DataFrame.

    Returns:
        str: The data strip with out the "---" row.
    """
    return re.sub(r"\n-{2,}\n", "\n", data)


def parse_median_intensity(file_name: Path) -> pd.DataFrame:
    return pd.read_csv(file_name)


def parse_imiss(file_name: Path):
    return pd.read_csv(file_name, delim_whitespace=True).rename({"IID": "Sample_ID"}, axis=1)


if __name__ == "__main__":
    if "snakemake" in locals():
        main(
            [Path(x) for x in snakemake.input.contamination],  # type: ignore # noqa
            Path(snakemake.input.median_idat_intensity),  # type: ignore # noqa
            Path(snakemake.input.imiss),  # type: ignore # noqa
            snakemake.params.intensity_threshold,  # type: ignore # noqa
            Path(snakemake.output[0]),  # type: ignore # noqa
        )
    else:
        app()

# Original code
# crDict3 = makeCallRateDict(input.imiss3)
# intensDict = {}
# with open(input.intens) as f:
#    head = f.readline()
#    line = f.readline()
#    while line != "":
#        (samp, chipId, intensity) = line.rstrip().split(",")
#        intensDict[samp] = float(intensity)
#        line = f.readline()
# with open(output[0], "w") as out:
#    out.write("ID,%Mix,LLK,LLK0\n")
#    for i in input.contam:
#        samp = os.path.basename(i).split(".")[0]
#        intens = intensDict[samp]
#        with open(i) as f:
#            head = f.readline()
#            while "%Mix" not in head and head != "":
#                head = f.readline()
#            if head == "":
#                print("strange file format: " + i)
#                sys.exit(1)
#            head = f.readline()
#            line = f.readline()
#            line_list = line.split()
#            line_list[0] = samp
#            if intens < params.intensThresh and not crDict3.get(samp):
#                line_list[1] = "NA"
#            out.write(",".join(line_list) + "\n")
