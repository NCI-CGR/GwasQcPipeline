import contextlib
import os
from hashlib import sha256
from math import isclose
from pathlib import Path
from textwrap import dedent
from typing import Optional, Union

from cgr_gwas_qc.config import config_to_yaml
from cgr_gwas_qc.models.config import (
    Config,
    ReferenceFiles,
    SoftwareParams,
    UserFiles,
    WorkflowParams,
)


@contextlib.contextmanager
def chdir(dirname: Path):
    curdir = Path().cwd()
    try:
        os.chdir(dirname)
        yield
    finally:
        os.chdir(curdir)


def file_hashes_equal(file1: Union[str, Path], file2: Union[str, Path]) -> bool:
    """Comapres two files using sha256 hashes.

    Takes two files and calculates their hashes returning True if they are
    the same and False otherwise.
    """
    return sha256(Path(file1).read_bytes()).digest() == sha256(Path(file2).read_bytes()).digest()


def file_rows_almost_equal(
    file1: Union[str, Path],
    file2: Union[str, Path],
    fuzzy_col: int,
    sep: str = "\t",
    header: bool = False,
) -> bool:
    """Compares two files row by row and makes sure they are almost equal"""
    results = []
    for row_idx, (r1, r2) in enumerate(
        zip(Path(file1).read_text().splitlines(), Path(file2).read_text().splitlines())
    ):
        if header and row_idx == 0:
            results.append(r1 == r2)
            continue

        for idx, (c1, c2) in enumerate(zip(r1.split(sep), r2.split(sep))):
            if idx != fuzzy_col:
                results.append(c1 == c2)
            else:
                results.append(isclose(float(c1), float(c2), rel_tol=1e-4))

    return all(results)


def make_test_config(
    current_dir: Path,
    data_type: str = "small",
    sample_sheet: str = "sample_sheet.csv",
    user_files: Optional[dict] = None,
) -> None:
    """Create a GwasQcPipeline test config.

    Args:
        current_dir: The working directory to save the `config.yml`.
        data_type: One of {"small", "real"}. Selects which kind of data to
          build the config for. Defaults to "small".
        sample_sheet: Name of the sample sheet. Defaults to "sample_sheet.csv".
        user_files: Dictionary of additional user files. If `None` then
          defines `gtc_pattern` and `idat_pattern` based on defaults.

    Raises:
        ValueError: if `data_type` is not "small" or "real".
    """

    if data_type == "small":
        project_name = "Small Test Project"
        reference_files = {
            "illumina_manifest_file": "small_manifest.bpm",
            "thousand_genome_vcf": "small_1KG.vcf.gz",
            "thousand_genome_tbi": "small_1KG.vcf.gz.tbi",
        }

        if user_files is None:
            user_files = {
                "gtc_pattern": "{Sample_ID}.gtc",
                "idat_pattern": {
                    "red": "{SentrixBarcode_A}_{SentrixPosition_A}_Red.idat",
                    "green": "{SentrixBarcode_A}_{SentrixPosition_A}_Grn.idat",
                },
            }

    elif data_type == "real":
        project_name = "Real Test Project"
        reference_files = {
            "illumina_manifest_file": "small_manifest.bpm",
            "thousand_genome_vcf": "small_1KG.vcf.gz",
            "thousand_genome_tbi": "small_1KG.vcf.gz.tbi",
        }

        if user_files is None:
            user_files = {
                "gtc_pattern": "{SentrixBarcode_A}_{SentrixPosition_A}.gtc",
                "idat_pattern": {
                    "red": "{SentrixBarcode_A}_{SentrixPosition_A}_Red.idat",
                    "green": "{SentrixBarcode_A}_{SentrixPosition_A}_Grn.idat",
                },
            }

    else:
        raise ValueError(
            f"`data_type` must be one of {{'small', 'real'}} and you entered: {data_type}"
        )

    # write out the config to `current_dir/config.yml`
    with chdir(current_dir):
        cfg = Config.construct(
            project_name=project_name,
            sample_sheet=sample_sheet,
            user_files=UserFiles.construct(**user_files),
            reference_files=ReferenceFiles.construct(**reference_files),
            software_params=SoftwareParams(),
            workflow_params=WorkflowParams(),
        )
        config_to_yaml(cfg)


def make_snakefile(working_dir: Path, contents: str):
    snakefile = working_dir / "Snakefile"
    snakefile.write_text(dedent(contents))
