"""Generate the internal sample QC report."""
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import typer

from cgr_gwas_qc import load_config
from cgr_gwas_qc.config import ConfigMgr
from cgr_gwas_qc.validators import check_file

app = typer.Typer(add_completion=False)


@app.command()
def main(
    imiss_start: Path = typer.Argument(..., help="Path to plink_start/samples.imiss"),
    imiss_cr1: Path = typer.Argument(..., help="Path to plink_filter_call_rate_1/samples.imiss"),
    imiss_cr2: Path = typer.Argument(..., help="Path to plink_filter_call_rate_2/samples.imiss"),
    sexcheck_cr1: Path = typer.Argument(
        ..., help="Path to plink_filter_call_rate_1/samples.sexcheck"
    ),
    ancestry: Path = typer.Argument(..., help="Path to ancestry/graf_ancestry_calls.txt"),
    known_concordance: Path = typer.Argument(..., help=""),
    unknown_concordance: Path = typer.Argument(..., help=""),
    # Optional inputs
    contam: Optional[Path] = typer.Option(
        None, help="Path to sample_filters/agg_contamination_test.csv"
    ),
    intensity: Optional[Path] = typer.Option(
        None, help="Path to sample_filters/agg_median_idat_intensity.csv"
    ),
    # Outputs
    all_qc: Path = typer.Argument("all_samples_qc.csv"),
    lims: str = typer.Argument("{out_name}_LimsUpload_{sample_sheet_date}.csv"),
):

    cfg = load_config()
    ss = _wrangle_sample_sheet(cfg.ss, cfg.config.workflow_params.expected_sex_col_name)
    Sample_IDs = ss.index

    df = pd.concat(
        [
            ss,
            _read_imiss_start(imiss_start, Sample_IDs),
            _read_imiss_cr1(imiss_cr1, Sample_IDs),
            _read_imiss_cr2(imiss_cr2, Sample_IDs),
            _read_sexcheck_cr1(sexcheck_cr1, ss.Expected_Sex),
            _read_ancestry(ancestry, Sample_IDs),
            _read_known_replicates(
                known_concordance, cfg.config.software_params.dup_concordance_cutoff, Sample_IDs
            ),
            _read_unknown_replicates(unknown_concordance, Sample_IDs),
            _read_contam(contam, cfg.config.software_params.contam_threshold, Sample_IDs),
            _read_intensity(intensity, Sample_IDs),
            _check_idats_files(cfg),
        ],
        axis=1,
    )

    qc_summary_flags = (
        "Low Call Rate",
        "Contaminated",
        "Sex Discordant",
        "Expected Replicate Discordance",
        "Unexpected Replicate",
    )
    df["Count_of_QC_Issue"] = df[qc_summary_flags].sum(axis=1)

    ################################################################################
    # Save Outputs
    ################################################################################
    _save_qc_table(df, all_qc)
    _save_lims_table(df, lims, cfg.config.sample_sheet)


def _wrangle_sample_sheet(df: pd.DataFrame, expected_sex_col_name: str) -> pd.DataFrame:
    """Identify expected sex column and count number samples per subject.

    Users can specify which column in the `sample_sheet` holds the expected
    sex information (`config.workflow_params.expected_sex_col_name`). While the QC
    report calls this column `Expected_Sex`. Here we rename
    `expected_sex_col_name` to "Expected_Sex".

    We also add the summary column with the number of `Sample_ID`s per
    `SR_Subject_ID`.

    Returns:
        pd.DataFrame: The full sample sheet with the following adjustments.
            - Sample_ID (pd.Index)
            - Expected_Sex (str): M/F based on the expected set column set in
              the config.
            - Count_of_SR_SubjectID (int): Number of `Sample_ID`s per
              `SR_SubjectID`.
    """
    return (
        df.rename({expected_sex_col_name: "Expected_Sex"}, axis=1)
        .merge(
            df.groupby("SR_Subject_ID", dropna=False).size().rename("Count_of_SR_SubjectID"),
            on="SR_Subject_ID",
        )  # Count the number of samples per subject ID.
        .set_index("Sample_ID")
    )


def _read_imiss_start(file_name: Path, Sample_IDs: pd.Index) -> pd.Series:
    """Read the starting call rates.

    Returns:
        pd.Series:
            - Sample_ID (pd.Index)
            - Call_Rate_Initial (float): The starting call rate calculated as
              `1 - F_MISS`.
    """
    return (
        pd.read_csv(file_name, delim_whitespace=True)  # FID IID MISS_PHENO N_MISS N_GENO F_MISS
        .rename({"IID": "Sample_ID"}, axis=1)
        .set_index("Sample_ID")
        .assign(Call_Rate_Initial=lambda x: 1 - x.F_MISS)
        .Call_Rate_Initial.squeeze()
        .reindex(Sample_IDs)
    )


def _read_imiss_cr1(file_name: Path, Sample_IDs: pd.Index) -> pd.DataFrame:
    """Read the call rates after level 1 filters.

    Returns:
        pd.DataFrame:
            - Sample_ID (pd.Index)
            - Call_Rate_1 (float): The call rate (`1 - F_MISS`) after level 1
              filters.
            - Call_Rate_1_filter (str): "Y" if call rate was less than the
              level 1 threshold else "N".
    """
    out_headers = ("Call_Rate_1", "Call_Rate_1_filter")
    return (
        pd.read_csv(file_name, delim_whitespace=True)  # FID IID MISS_PHENO N_MISS N_GENO F_MISS
        .rename({"IID": "Sample_ID"}, axis=1)
        .set_index("Sample_ID")
        .assign(Call_Rate_1=lambda x: 1 - x.F_MISS)
        .reindex(Sample_IDs)
        .assign(Call_Rate_1_filter=lambda x: np.where(x.Call_Rate_1.isnull(), "Y", "N"))
        .loc[:, out_headers]
    )


def _read_imiss_cr2(file_name: Path, Sample_IDs: pd.Index) -> pd.DataFrame:
    """Read the call rates after level 2 filters.

    Returns:
        pd.DataFrame:
            - Sample_ID (pd.Index)
            - Call_Rate_2 (float): The call rate (`1 - F_MISS`) after level 2
              filters.
            - Call_Rate_2_filter (str): "Y" if call rate was less than the
              level 2 threshold else "N".
    """
    out_headers = ("Call_Rate_2", "Call_Rate_2_filter", "Low Call Rate")
    return (
        pd.read_csv(file_name, delim_whitespace=True)  # FID IID MISS_PHENO N_MISS N_GENO F_MISS
        .rename({"IID": "Sample_ID"}, axis=1)
        .set_index("Sample_ID")
        .assign(Call_Rate_2=lambda x: 1 - x.F_MISS)
        .reindex(Sample_IDs)
        .assign(Call_Rate_2_filter=lambda x: np.where(x.Call_Rate_2.isnull(), "Y", "N"))
        .assign(**{"Low Call Rate": lambda x: x.Call_Rate_2_filter == "Y"})
        .loc[:, out_headers]
    )


def _read_sexcheck_cr1(file_name: Path, Expected_Sex: pd.Series) -> pd.DataFrame:
    """Read sex predictions and summarize.

    Read PLINK sex prediction file. Convert the `Predicted_Sex` indicator
    variable to M/F designations. Compare predicted results with the expected
    sexes and create a summary column `SexMatch` if predicted/expected sex
    calls match. Then flag samples as `Sex Discordant` if sex was predicted
    to be different than expected.

    Returns:
        pd.DataFrame:
            - Sample_ID (pd.Index)
            - ChrX_Inbreed_estimate (float): PLINK's inbreeding coefficient
              from sexcheck.
            - Predicted_Sex (str): M/F/U based on PLINK sex predictions.
            - SexMatch (str): Y if expected and predicted are equal. N if they
              are different. U if prediction was U.
            - Sex Discordant (bool): True if SexMatch == "N"
    """
    df = (
        pd.read_csv(file_name, delim_whitespace=True)
        .rename({"IID": "Sample_ID", "F": "ChrX_Inbreed_estimate"}, axis=1)
        .set_index("Sample_ID")
        .assign(Predicted_Sex=lambda x: x.SNPSEX.map({1: "M", 2: "F"}))
        .reindex(Expected_Sex.index)
        .fillna({"Predicted_Sex": "U"})
        .loc[:, ["ChrX_Inbreed_estimate", "Predicted_Sex"]]
    )

    # TODO: Decide if we want to keep this logic. See
    # http://10.133.130.114/jfear/GwasQcPipeline/issues/35
    df.loc[df.ChrX_Inbreed_estimate < 0.5, "Predicted_Sex"] = "F"
    df.loc[df.ChrX_Inbreed_estimate >= 0.5, "Predicted_Sex"] = "M"

    df["SexMatch"] = np.where(Expected_Sex == df.Predicted_Sex, "Y", "N")
    df.loc[df.Predicted_Sex == "U", "SexMatch"] = "U"  # If we could not predict sex then label as U

    df["Sex Discordant"] = df.SexMatch.replace({"N": True, "Y": False, "U": np.nan})

    return df


def _read_ancestry(file_name: Path, Sample_IDs: pd.Index) -> pd.DataFrame:
    """Read GRAF's ancestry calls.

    Returns:
        pd.DataFrame:
            - Sample_ID (pd.Index)
            - AFR (float): The proportion African ancestry (`P_f(%) / 100`).
            - EUR (float): The proportion African ancestry (`P_e(%) / 100`).
            - ASN (float): The proportion African ancestry (`P_a(%) / 100`).
            - Ancestry (str): The ancestry label assigned by GRAF. See tables
              2 and 3 in their manuscript_ for the exact criteria.

    .. _manuscript: https://pubmed.ncbi.nlm.nih.gov/31151998/

    """
    return (
        pd.read_csv(file_name, sep="\t")
        .rename({"Subject": "Sample_ID", "Computed population": "Ancestry"})
        .assign(AFR=lambda x: x["P_f (%)"] / 100)
        .assign(EUR=lambda x: x["P_e (%)"] / 100)
        .assign(ASN=lambda x: x["P_a (%)"] / 100)
        .set_index("Sample_ID")
        .loc[:, ("AFR", "EUR", "ASN", "Ancestry")]
        .reindex(Sample_IDs)
    )


def _read_known_replicates(
    file_name: Path, dup_concordance_cutoff: float, Sample_IDs: pd.Index
) -> pd.Series:  # noqa
    """Create a flag of known replicates that show low concordance.

    Given a set of samples that are known to be from the same Subject. Flag
    samples that show low concordance with one or more replicates.

    Returns:
        pd.Series:
            - Sample_ID (pd.Index)
            - Expected Replicate Discordance (bool): True if replicates show
              a concordance below the supplied threshold. Otherwise False.
   """

    discord_Sample_IDs = (
        pd.read_csv(file_name)  # Subject_ID Sample_ID1 Sample_ID2 Concordance PI_HAT
        .query("Concordance.notna() & Concordance < @dup_concordance_cutoff")
        .loc[:, ("Sample_ID1", "Sample_ID2")]
        .melt()
        .value.unique()
    )  # A set of Sample_IDs that were replicates were not concordant.

    sr = pd.Series(False, index=Sample_IDs).rename("Expected Replicate Discordance")
    sr[sr.index.isin(discord_Sample_IDs)] = True

    return sr


def _read_unknown_replicates(file_name: Path, Sample_IDs: pd.Index) -> pd.Series:
    """Create a flag of for samples that are likely to be unknown replicates.

    Given a set of samples that are not known to be from the same Subject but
    show high concordance.

    Returns:
        pd.Series:
            - Sample_ID (pd.Index)
            - Expected Replicate Discordance (bool): True if replicates show
              a concordance below the supplied threshold. Otherwise False.
   """

    cord_Sample_IDs = (
        pd.read_csv(file_name)  # Subject_ID1 Subject_ID2 Sample_ID1 Sample_ID2 Concordance PI_HAT
        .loc[:, ("Sample_ID1", "Sample_ID2")]
        .melt()
        .value.unique()
    )  # A set of Sample_IDs that look like a replicate with another sample from a different subject.

    sr = pd.Series(False, index=Sample_IDs).rename("Unexpected Replicate")
    sr[sr.index.isin(cord_Sample_IDs)] = True

    return sr

    pass


def _read_contam(
    file_name: Optional[Path], contam_threshold: float, Sample_IDs: pd.Index
) -> pd.DataFrame:
    """Parse verifyIDintensity contamination information.

    Returns:
        pd.DataFrame:
            - Sample_ID (pd.index)
            - Contamination_Rate (float): The contamination rate as estimated
              by verifyIDintensity.
            - Contaminated (bool): True if the contamination rate is greater
              than the supplied threshold.
    """

    if file_name is None:
        return pd.DataFrame(index=Sample_IDs, columns=["Contamination_Rate", "Contaminated"])

    df = (
        pd.read_csv(file_name)
        .rename({"ID": "Sample_ID", "%Mix": "Contamination_Rate"})
        .set_index("Sample_ID")
    )

    df["Contaminated"] = df.Contamination_Rate > contam_threshold
    df.loc[df.Contamination_Rate.isna(), "Contaminated"] = np.nan

    return df.reindex(Sample_IDs)[("Contamination_Rate", "Contaminated")]


def _read_intensity(file_name: Optional[Path], Sample_IDs: pd.Index) -> pd.Series:
    """Parse the median Idat intensity table.

    Returns:
        pd.Series:
            - Sample_ID (pd.Index)
            - IdatIntensity (float): The median Idat intensity.
    """
    if file_name is None:
        return pd.Series(index=Sample_IDs).rename("IdatIntensity")

    return (
        pd.read_csv(file_name)
        .rename({"SampId": "Sample_ID", "MedianIntensity": "IdatIntensity"})
        .set_index("Sample_ID")
        .reindex(Sample_IDs)
        .IdatIntensity
    )


def _check_idats_files(cfg: ConfigMgr) -> pd.Series:
    """Check that red and green IDAT files exist.

    Args:
        df: A sample table with at least `Sampel_ID` and columns needed to fill wildcards.
        red: A wildcard pattern for red files. Wildcards must be in `df`.
        green: A wildcard pattern for green files. Wildcards must be in `df`.

    Returns:
        pd.Series:
            - Sample_ID (pd.Index)
            - IdatsInProjectDir (bool): True if both the red and green Idat files existed
    """
    if cfg.config.user_files.idat_pattern is None:
        # No Idat path specified in config, return all NaN.
        return pd.Series(index=cfg.ss.Sample_ID).rename("IdatsInProjectDir")

    results = []
    for Sample_ID, red, green in zip(
        cfg.ss.Sample_ID,
        cfg.expand(cfg.config.user_files.idat_pattern.red),
        cfg.expand(cfg.config.user_files.idat_pattern.green),
    ):
        try:
            check_file(Path(red))
            check_file(Path(green))
            results.append((Sample_ID, "YES"))
        except (FileNotFoundError, PermissionError):
            results.append((Sample_ID, "NO"))

    return pd.Series(dict(results)).rename_axis("Sample_ID").rename("IdatsInProjectDir")


def _save_qc_table(df: pd.DataFrame, file_name: Path) -> None:
    header_order = [
        "SR_Subject_ID",
        "Count_of_SR_SubjectID",
        "SR",
        "Current_Subject_Status",
        "Subject_Notes",
        "LIMS_Individual_ID",
        "Project",
        "Sample_ID",
        "Project-Sample ID",
        "LIMSSample_ID",
        "Sample_Status",
        "IdatsInProjectDir",
        "IdatIntensity",
        "Expected_Sex",
        "Predicted_Sex",
        "SexMatch",
        "ChrX_Inbreed_estimate",
        "AFR",
        "EUR",
        "ASN",
        "Ancestry",
        "Contamination_Rate",
        "Call_Rate_Initial",
        "Call_Rate_1_filter",
        "Call_Rate_1",
        "Call_Rate_2_filter",
        "Call_Rate_2",
        "Low Call Rate",
        "Contaminated",
        "Sex Discordant",
        "Expected Replicate Discordance",
        "Unexpected Replicate",
        "Count_of_QC_Issue",
        "Identifiler_Needed",
    ]

    df.reset_index().loc[:, header_order].to_csv(file_name, index=False)


def _save_lims_table(df: pd.DataFrame, file_name_pattern: str, sample_sheet_path: Path) -> None:
    if "AnalysisManifest" in sample_sheet_path.name:
        # example file_name_pattern: "{out_name}_LimsUpload_{sample_sheet_date}.csv"
        # example sample_sheet name: "SR0001-001_1_AnalysisManifest_999999999999_20200101.csv"
        out_name, sample_sheet_date = sample_sheet_path.stem.split("_AnalysisManifest_")

        # example file_name: "SR0001-001_1_LimsUpload_999999999999_20200101.csv"
        file_name = file_name_pattern.format(out_name=out_name, sample_sheet_date=sample_sheet_date)
    else:
        file_name = f"{sample_sheet_path.stem}_LimsUpload.csv"

    header_order = [
        "SR_Subject_ID",
        "LIMS_Individual_ID",
        "Sample ID",
        "Project-Sample ID",
        "Call Rate",
        "Low Call Rate",
        "Contaminated",
        "Sex Discordant",
        "Expected Replicate Discordance",
        "Unexpected Replicate",
    ]

    df.reset_index().loc[:, header_order].to_csv(file_name, index=False)


if __name__ == "__main__":
    if "snakemake" in locals():
        main(**snakemake.input, **snakemake.output)  # type: ignore # noqa
    else:
        app()
