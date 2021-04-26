from dataclasses import dataclass
from pathlib import Path
from string import ascii_lowercase
from typing import List

import pandas as pd


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
            SexVerification.construct(subject_qc, chrx_inbreeding_png),
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

    @classmethod
    def construct(cls, subject_qc: pd.DataFrame, png: Path) -> "SexVerification":
        return cls(
            subject_qc.is_sex_discordant.sum(),
            (~subject_qc.is_sex_discordant).sum(),
            png.resolve().as_posix(),
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
