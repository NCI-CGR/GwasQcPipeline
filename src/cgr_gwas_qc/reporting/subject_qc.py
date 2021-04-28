from dataclasses import dataclass
from pathlib import Path
from string import ascii_lowercase
from typing import List, Optional

import pandas as pd

from cgr_gwas_qc.reporting.constants import REPORT_NAME_MAPPER


@dataclass
class SubjectQC:
    unexpected_replicates: "UnExpectedReplicates"
    sex_verification: "SexVerification"
    relatives: "Relatedness"
    autosomal: "Autosomal"
    pca: "Pca"
    hwe: "Hwe"

    @classmethod
    def construct(
        cls,
        ss: pd.DataFrame,
        subject_qc: pd.DataFrame,
        unexpected_replicates: Path,
        chrx_inbreeding_png: Path,
        population_qc: pd.DataFrame,
        autosomal_heterozygosity_png_dir: Path,
        pca_png_dir: Path,
        hwe_png_dir: Path,
    ) -> "SubjectQC":
        return cls(
            UnExpectedReplicates.construct(subject_qc, unexpected_replicates),
            SexVerification.construct(ss, subject_qc, chrx_inbreeding_png),
            Relatedness.construct(population_qc),
            Autosomal.construct(population_qc, autosomal_heterozygosity_png_dir),
            Pca.construct(pca_png_dir),
            Hwe.construct(hwe_png_dir),
        )


@dataclass
class UnExpectedReplicates:
    num_unexpected_replicates: int
    min_concordance: float
    mean_concordance: float

    @classmethod
    def construct(
        cls, subject_qc: pd.DataFrame, replicates: pd.DataFrame
    ) -> "UnExpectedReplicates":
        return cls(
            subject_qc.is_unexpected_replicate.sum(),
            replicates.PLINK_concordance.min(),
            replicates.PLINK_concordance.mean(),
        )


@dataclass
class SexVerification:
    num_sex_discordant: int
    num_remaining: int
    png: str
    table: Optional[str]

    @classmethod
    def construct(cls, ss: pd.DataFrame, subject_qc: pd.DataFrame, png: Path) -> "SexVerification":
        return cls(
            subject_qc.is_sex_discordant.sum(),
            (~subject_qc.is_sex_discordant).sum(),
            png.resolve().as_posix(),
            cls.build_table(ss, subject_qc),
        )

    @staticmethod
    def build_table(ss: pd.DataFrame, subject_qc: pd.DataFrame) -> Optional[str]:
        cols = [
            "Sample_ID",
            "Group_By_Subject_ID",
            "Sample_Name",
            "Identifiler_Sex",
            "Project",
            "expected_sex",
            "predicted_sex",
            "XY characters",
        ]
        df = subject_qc.query("is_sex_discordant")

        if df.shape[0] == 0:
            return None

        return (
            df.merge(ss, how="left", on="Sample_ID", suffixes=["", "_DROP"])
            .reindex(cols, axis=1)
            .rename(REPORT_NAME_MAPPER, axis=1)
            .set_index("Sample_ID")
            .fillna({"XY characters": "Not Analyzed"})
            .to_markdown()
        )


@dataclass
class Relatedness:
    num_related_subjects: int
    num_qc_families: int

    @classmethod
    def construct(cls, population_qc: pd.DataFrame) -> "Relatedness":
        return cls(
            population_qc.relatives.dropna().shape[0],
            population_qc.QC_Family_ID.dropna().unique().shape[0],
        )


@dataclass
class FigurePanel:
    letter: str
    population: str
    png: str

    @classmethod
    def from_png_dir(cls, dirname: Path) -> List["FigurePanel"]:
        panels = []
        for letter, filename in zip(
            ascii_lowercase, sorted(dirname.glob("*.png"), key=lambda x: x.stem)
        ):
            population = filename.stem
            panels.append(cls(letter, population, filename.resolve().as_posix()))

        return panels


@dataclass
class Autosomal:
    num_populations_analyzed: int
    num_subjects_excluded: int
    panels: List[FigurePanel]

    @classmethod
    def construct(cls, population_qc: pd.DataFrame, png_dir: Path) -> "Autosomal":
        return cls(
            population_qc.population.unique().shape[0],
            population_qc.is_extreme_autosomal_heterozygosity.sum(),
            panels=FigurePanel.from_png_dir(png_dir),
        )


@dataclass
class Pca:
    panels: List[FigurePanel]

    @classmethod
    def construct(cls, png_dir: Path) -> "Pca":
        return cls(panels=FigurePanel.from_png_dir(png_dir))


@dataclass
class Hwe:
    panels: List[FigurePanel]

    @classmethod
    def construct(cls, png_dir: Path) -> "Hwe":
        return cls(panels=FigurePanel.from_png_dir(png_dir))
