from dataclasses import dataclass
from pathlib import Path
from string import ascii_lowercase
from typing import List

import pandas as pd


@dataclass
class SubjectQC:
    relatives: "Relatedness"
    autosomal: "Autosomal"
    pca: "Pca"
    hwe: "Hwe"

    @classmethod
    def construct(
        cls,
        population_qc: pd.DataFrame,
        autosomal_heterozygosity_png_dir: Path,
        pca_png_dir: Path,
        hwe_png_dir: Path,
    ) -> "SubjectQC":
        return cls(
            relatives=Relatedness.construct(population_qc),
            autosomal=Autosomal.construct(population_qc, autosomal_heterozygosity_png_dir),
            pca=Pca.construct(pca_png_dir),
            hwe=Hwe.construct(hwe_png_dir),
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
