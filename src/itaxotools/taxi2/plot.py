from enum import Enum
from pathlib import Path
from typing import NamedTuple

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

matplotlib.use("agg")
matplotlib.use("svg")
matplotlib.use("pdf")


class ComparisonType(Enum):
    Unknown = "no info"
    IntraSpecies = "intra-species"
    InterSpecies = "inter-species"
    IntraGenus = "intra-genus"
    InterGenus = "inter-genus"

    def __init__(self, label):
        self.index = len(type(self).__members__)
        self.label = label

    def __lt__(self, other):
        return self.index < other.index


class HistogramPoint(NamedTuple):
    type: str
    value: float


class HistogramPlotter:
    def __init__(
        self, formats: list[str] = None, palette=None, binwidth=0.05, binfactor=1.0
    ):
        self.formats = formats or ["png", "svg", "pdf"]
        self.palette = palette or sns.color_palette()
        self.binwidth = binwidth
        self.binfactor = binfactor

        self.metrics: dict[str, HistogramPoint] = dict()

    def add(self, metric: str, value: float, type: ComparisonType):
        if metric not in self.metrics:
            self.metrics[metric] = list()
        point = HistogramPoint(type.label, value)
        self.metrics[metric].append(point)

    def plot(self, output_path: Path):
        matplotlib.rcParams["pdf.fonttype"] = 42
        matplotlib.rcParams["ps.fonttype"] = 42
        matplotlib.rc("font", **{"family": "sans-serif"})

        for metric in self.metrics:
            path = output_path / metric
            path.mkdir(exist_ok=True)
            df = pd.DataFrame(self.metrics[metric], columns=HistogramPoint._fields)
            types = df["type"].unique()
            palette, order = self.palette_from_types(types)

            self.plot_layered(
                metric, df, palette, order, path / f"{metric}_layered_hist"
            )
            self.plot_histogram(
                metric, df, "stack", palette, order, path / f"{metric}_stacked_hist"
            )
            self.plot_histogram(
                metric, df, "dodge", palette, order, path / f"{metric}_dodge_hist"
            )

            has_species_info = (
                ComparisonType.IntraSpecies.label in types
                or ComparisonType.InterSpecies.label in types
            )
            has_genus_info = (
                ComparisonType.IntraGenus.label in types
                or ComparisonType.InterGenus.label in types
            )

            if has_species_info and has_genus_info:
                # species info only

                translation = {
                    ComparisonType.InterGenus: ComparisonType.InterSpecies,
                    ComparisonType.IntraGenus: ComparisonType.Unknown,
                }

                df_species = df.copy()
                for old, new in translation.items():
                    df_species["type"] = df_species["type"].str.replace(
                        old.label, new.label, regex=False
                    )

                types = df_species["type"].unique()
                palette, order = self.palette_from_types(types)

                path_species = path / "species_only"
                path_species.mkdir(exist_ok=True)

                self.plot_layered(
                    metric,
                    df_species,
                    palette,
                    order,
                    path_species / f"{metric}_layered_hist_species_only",
                )
                self.plot_histogram(
                    metric,
                    df_species,
                    "stack",
                    palette,
                    order,
                    path_species / f"{metric}_stacked_hist_species_only",
                )
                self.plot_histogram(
                    metric,
                    df_species,
                    "dodge",
                    palette,
                    order,
                    path_species / f"{metric}_dodge_hist_species_only",
                )

                # genus info only

                translation = {
                    ComparisonType.InterSpecies: ComparisonType.IntraGenus,
                    ComparisonType.IntraSpecies: ComparisonType.IntraGenus,
                }

                df_genus = df.copy()
                for old, new in translation.items():
                    df_genus["type"] = df_genus["type"].str.replace(
                        old.label, new.label, regex=False
                    )

                types = df_genus["type"].unique()
                palette, order = self.palette_from_types(types)

                path_genus = path / "genus_only"
                path_genus.mkdir(exist_ok=True)

                self.plot_layered(
                    metric,
                    df_genus,
                    palette,
                    order,
                    path_genus / f"{metric}_layered_hist_genus_only",
                )
                self.plot_histogram(
                    metric,
                    df_genus,
                    "stack",
                    palette,
                    order,
                    path_genus / f"{metric}_stacked_hist_genus_only",
                )
                self.plot_histogram(
                    metric,
                    df_genus,
                    "dodge",
                    palette,
                    order,
                    path_genus / f"{metric}_dodge_hist_genus_only",
                )

    def plot_layered(
        self,
        metric: str,
        df: pd.DataFrame,
        palette: list[tuple],
        order: list[str],
        path: Path,
    ):
        g = sns.FacetGrid(
            df,
            row="type",
            hue="type",
            palette=palette,
            hue_order=order,
            height=1.5,
            aspect=4,
        )
        g.map_dataframe(
            sns.histplot,
            x="value",
            binwidth=self.binwidth * self.binfactor,
            binrange=(0.0, self.binfactor),
        )
        g.set_xlabels(f"{metric} distance")
        g.set_ylabels("Count")
        for format in self.formats:
            g.savefig(path.with_suffix(f".{format}"), transparent=True)
        plt.close(g.fig)

    def plot_histogram(
        self,
        metric: str,
        df: pd.DataFrame,
        multiple: str,
        palette: list[tuple],
        order: list[str],
        path: Path,
    ):
        fig, ax = plt.subplots()
        sns.histplot(
            df,
            x="value",
            hue="type",
            multiple=multiple,
            binwidth=self.binwidth * self.binfactor,
            binrange=(0.0, self.binfactor),
            palette=palette,
            hue_order=order,
            ax=ax,
        )
        sns.despine()

        ax.set_xlabel(f"{metric} distance")
        ax.set_ylabel("Count")

        for format in self.formats:
            fig.savefig(path.with_suffix(f".{format}"), transparent=True)
        plt.close(fig)

    def palette_from_types(self, types: list[str]) -> tuple[list[tuple], list[str]]:
        """Make sure each type has consistent color among runs"""
        types = sorted([ComparisonType(type) for type in types])
        palette = [self.palette[type.index] for type in types]
        order = [type.label for type in types]
        return palette, order
