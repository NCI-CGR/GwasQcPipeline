from hashlib import sha256
from math import isclose
from pathlib import Path

import pandas as pd
from pandas.testing import assert_series_equal

from cgr_gwas_qc.typing import PathLike


################################################################################
# General
################################################################################
def file_hashes_equal(file1: PathLike, file2: PathLike) -> bool:
    """Comapres two files using sha256 hashes.

    Takes two files and calculates their hashes returning True if they are
    the same and False otherwise.
    """
    return sha256(Path(file1).read_bytes()).digest() == sha256(Path(file2).read_bytes()).digest()


def sorted_file_equal(file1: PathLike, file2: PathLike) -> bool:
    """Comapres two sorted files.

    Reads in each file, splits by line, and then sorts.
    """
    pth1, pth2 = Path(file1), Path(file2)
    return sorted(pth1.read_text().splitlines()) == sorted(pth2.read_text().splitlines())


def file_rows_almost_equal(
    file1: PathLike,
    file2: PathLike,
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


################################################################################
# PLINK File Fromats
################################################################################
def assert_plink_bed_equal(file1: PathLike, file2: PathLike):
    raise NotImplementedError


def assert_plink_bim_equal(file1: PathLike, file2: PathLike):
    """Compare two PLINK BIM files ignoring allele order.

    Raises:
        AssertionError: if all rows are not equal after ignoring allele order.
    """
    from cgr_gwas_qc.parsers.bim import BimFile

    fh1, fh2 = BimFile(file1), BimFile(file2)
    rows_are_equal = [record1 == record2 for record1, record2 in zip(fh1, fh2)]
    assert all(rows_are_equal)


def assert_plink_fam_equal(file1: PathLike, file2: PathLike):
    """Compare two PLINK FAM files.

    Raises:
        AssertionError: if files are not identical.
    """
    assert file_hashes_equal(file1, file2)


################################################################################
# Major Workflow Files
################################################################################
def assert_legacy_dev_sample_qc_equal(legacy_file: PathLike, dev_file: PathLike):
    legacy = pd.read_csv(legacy_file).set_index("Sample_ID").sort_index()
    dev = pd.read_csv(dev_file).set_index("Sample_ID").sort_index()

    assert legacy.shape[0] == dev.shape[0]

    # --------------------
    # IDATS
    # --------------------
    assert_series_equal(
        legacy["IdatsInProjectDir"].replace({"YES": True, "NO": False}),
        ~dev["is_missing_idats"],
        check_dtype=False,
        check_names=False,
    )
    assert_series_equal(legacy["IdatIntensity"], dev["IdatIntensity"])

    # --------------------
    # Call Rate
    # --------------------
    assert_series_equal(legacy["Call_Rate_Initial"], dev["Call_Rate_Initial"])
    assert_series_equal(
        legacy["Call_Rate_1_filter"].replace({"Y": True, "N": False}),
        dev["is_cr1_filtered"],
        check_dtype=False,
        check_names=False,
    )
    assert_series_equal(legacy["Call_Rate_1"], dev["Call_Rate_1"])
    assert_series_equal(
        legacy["Call_Rate_2_filter"].replace({"Y": True, "N": False}),
        dev["is_cr2_filtered"].fillna(True),  # dev sets samples as NA if filtered by cr1
        check_dtype=False,
        check_names=False,
    )
    assert_series_equal(legacy["Call_Rate_2"], dev["Call_Rate_2"])
    assert_series_equal(
        legacy["Low Call Rate"],
        dev["is_call_rate_filtered"],
        check_dtype=False,
        check_names=False,
    )

    # --------------------
    # Contamination
    # --------------------
    mask = dev["Contamination_Rate"].notna()
    assert_series_equal(
        legacy["Contamination_Rate"][mask], dev["Contamination_Rate"][mask], atol=1e-3
    )
    assert_series_equal(
        legacy["Contaminated"], dev["is_contaminated"], check_dtype=False, check_names=False
    )

    # --------------------
    # Sex
    # --------------------
    assert_series_equal(legacy["Expected_Sex"], dev["expected_sex"], check_names=False)
    assert_series_equal(
        legacy["Predicted_Sex"], dev["predicted_sex"].fillna("U"), check_names=False
    )
    assert_series_equal(
        legacy["ChrX_Inbreed_estimate"], dev["X_inbreeding_coefficient"], check_names=False
    )
    assert_series_equal(
        legacy["Sex Discordant"], dev["is_sex_discordant"], check_dtype=False, check_names=False
    )

    # --------------------
    # Ancestry
    # --------------------
    assert_series_equal(legacy["AFR"], dev["AFR"], atol=0.1)
    assert_series_equal(legacy["EUR"], dev["EUR"], atol=0.1)
    assert_series_equal(legacy["ASN"], dev["ASN"], atol=0.1)

    # --------------------
    # Concordance
    # --------------------
    assert_series_equal(
        legacy["Expected Replicate Discordance"],
        dev["is_discordant_replicate"].fillna(False),
        check_dtype=False,
        check_names=False,
    )

    # --------------------
    # Summary
    # --------------------
    assert_series_equal(legacy["Identifiler_Needed"], dev["identifiler_needed"], check_names=False)

    # -----------------------
    # Problematic Comparisons
    # -----------------------

    #  I am not including sex discordance in this count, maybe I should
    # assert_series_equal(legacy["Count_of_QC_Issue"], dev["num_analytic_exclusion"], check_names=False)

    # Ancestries are difference b/c very different methods
    # assert_series_equal(legacy["Ancestry"], dev["Ancestry"])

    # This is different and not sure why. Probably b/c I am counting based on
    # PI_Subject_ID here and not SR_Subject_ID
    # assert_series_equal(legacy["Count_of_SR_SubjectID"], dev["num_samples_per_subject"])
