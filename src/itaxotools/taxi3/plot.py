from typing import NamedTuple
from pathlib import Path

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc


class HistogramPoint(NamedTuple):
    _types = [
        'intra-species',
        'inter-species',
        'intra-genus',
        'inter-genus',
        '-',
    ]

    value: float
    type: str


class HistogramPlotter:
    def __init__(self, formats: list[str] = None, palette = None, binwidth = 0.05, binfactor = 1.0):
        self.formats = formats or ['png']
        self.palette = palette or sns.color_palette()
        self.binwidth = binwidth
        self.binfactor = binfactor

        self.metrics: dict[str, HistogramPoint] = dict()

    def add(self, metric: str, value: float, type: str):
        if not metric in self.metrics:
            self.metrics[metric] = list()
        point = HistogramPoint(value, type)
        self.metrics[metric].append(point)

    def plot(self, output_path: Path):
        matplotlib.use('agg')
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42
        rc('font', **{'family': 'sans-serif'})

        for metric in self.metrics:
            path = output_path / metric
            path.mkdir(exist_ok=True)
            df = pd.DataFrame(self.metrics[metric], columns=HistogramPoint._fields)
            types = list(df['type'].unique())
            palette = self.palette_from_types(types)

            self.plot_layered(metric, df, palette, path / f'{metric}_layered_hist')
            self.plot_histogram(metric, df, 'stack', palette, path / f'{metric}_stacked_hist')
            self.plot_histogram(metric, df, 'dodge', palette, path / f'{metric}_dodge_hist')

            if not all((
                any(type in types for type in ['inter-species', 'intra-species']),
                any(type in types for type in ['inter-genus', 'intra-genus']),
            )):
                continue

            path = path / 'no_genus'
            path.mkdir(exist_ok=True)
            df['type'] = df['type'].str.replace('inter-genus', 'inter-species')
            if 'inter-genus' in types:
                types.remove('inter-genus')
            palette = self.palette_from_types(types)
            self.plot_layered(metric, df, palette, path / f'{metric}_layered_hist_no_genus')
            self.plot_histogram(metric, df, 'stack', palette, path / f'{metric}_stacked_hist_no_genus')
            self.plot_histogram(metric, df, 'dodge', palette, path / f'{metric}_dodge_hist_no_genus')

    def plot_layered(self, metric: str, df: pd.DataFrame, palette: list[tuple], path: Path):
        g = sns.FacetGrid(df, row='type', hue='type', palette=palette, height=1.5, aspect=4)
        g.map_dataframe(sns.histplot, x='value', binwidth=self.binwidth * self.binfactor, binrange=(0.0, self.binfactor))
        g.set_xlabels(f'{metric} distance')
        g.set_ylabels('Count')
        for format in self.formats:
            g.savefig(path.with_suffix(f'.{format}'), transparent=True)
        plt.close(g.fig)

    def plot_histogram(self, metric: str, df: pd.DataFrame, multiple: str, palette: list[tuple], path: Path):
        fig, ax = plt.subplots()
        sns.histplot(df, x='value', hue='type', multiple=multiple, binwidth=self.binwidth * self.binfactor, binrange=(0.0, self.binfactor), palette=palette, ax=ax)
        sns.despine()

        ax.set_xlabel(f'{metric} distance')
        ax.set_ylabel('Count')

        for format in self.formats:
            fig.savefig(path.with_suffix(f'.{format}'), transparent=True)
        plt.close(fig)

    def palette_from_types(self, types: list[str]) -> list[tuple]:
        """Make sure each type has consistent color among runs"""
        indices = [HistogramPoint._types.index(type) for type in types]
        palette = self.palette
        return [palette[index] for index in indices]
