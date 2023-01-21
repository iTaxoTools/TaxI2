from typing import NamedTuple
from pathlib import Path

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc


class HistogramPoint(NamedTuple):
    value: float
    type: str


class HistogramPlotter:
    def __init__(self, formats: list[str] = None, binwindth = 0.05):
        self.groups: dict[str, HistogramPoint] = dict()
        self.formats = formats or ['png']
        self.binwindth = binwindth

    def add(self, group: str, value: float, type: str):
        if not group in self.groups:
            self.groups[group] = list()
        point = HistogramPoint(value, type)
        self.groups[group].append(point)

    def plot(self, path: Path):
        matplotlib.use("agg")
        matplotlib.rcParams["pdf.fonttype"] = 42
        matplotlib.rcParams["ps.fonttype"] = 42
        rc("font", **{"family": "sans-serif"})

        self.plot_layered(path / 'layered')

    def plot_layered(self, path: Path):
        path.mkdir(parents=True, exist_ok=True)
        for group, points in self.groups.items():
            df = pd.DataFrame(self.groups[group], columns=HistogramPoint._fields)
            g = sns.FacetGrid(df, row='type', height=1.5, aspect=4)
            g.map_dataframe(sns.histplot, x='value', binwidth=self.binwindth, binrange=(0.0, 1.0))
            g.set_axis_labels('distance', 'Count')
            g.set_xlabels(group)
            for format in self.formats:
                g.savefig(path / f'{group}_layered_hist.{format}', transparent=True)
            plt.close(g.fig)
