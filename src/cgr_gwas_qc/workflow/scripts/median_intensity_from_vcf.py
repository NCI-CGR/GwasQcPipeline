#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from pathlib import Path
from statistics import median

import numpy as np
import pandas as pd
import typer
from pysam import VariantFile

app = typer.Typer(add_completion=False)


@app.command()
def main(
    vcf_file: Path = typer.Argument(
        ...,
        help="Path to a multisample BCF/VCF file with needed scores.",
        exists=True,
        readable=True,
    ),
    sample_id: str = typer.Argument(..., help="Sample_ID. Should match as in VCF"),
    outfile: Path = typer.Argument(
        ...,
        help="Path to output file to write per sample median intensity.",
        file_okay=True,
        writable=True,
    ),
) -> None:
    fields_to_fetch = ["X", "Y"]
    vcf = VariantFile(vcf_file)

    # loops over all records in the VCF, fetches the 'X' and 'Y' intensity values, sums them and stores them in the list
    # of 'sample_values'.z
    sample_values = []
    for rec in vcf:
        sample_values.append(np.sum([rec.samples[sample_id].get(key) for key in fields_to_fetch]))
    vcf.close()

    median_intensity = median(sample_values)
    df = pd.DataFrame({"Sample_ID": [sample_id], "median_intensity": [median_intensity]})
    df.to_csv(outfile, index=False)


if __name__ == "__main__" and "snakemake" in locals():
    main(
        **{k: v for k, v in snakemake.input.items() if not k.startswith("_")},  # type: ignore # noqa
        **{k: v for k, v in snakemake.params.items()},  # type: ignore # noqa
        outfile=snakemake.output[0],  # type: ignore # noqa
    )
