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
    def __init__(self, formats: list[str] = None, binwindth = 0.05, palette = None):
        self.metrics: dict[str, HistogramPoint] = dict()
        self.formats = formats or ['png']
        self.binwindth = binwindth
        self.palette = palette

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
            self.plot_layered(metric, df, path / f'{metric}_layered_hist')
            self.plot_histogram(metric, df, 'stack', path / f'{metric}_stacked_hist')
            self.plot_histogram(metric, df, 'dodge', path / f'{metric}_dodge_hist')
            df['type'] = df['type'].str.replace('inter-genus', 'inter-species')
            self.plot_layered(metric, df, path / f'{metric}_layered_hist_no_genus')
            self.plot_histogram(metric, df, 'stack', path / f'{metric}_stacked_hist_no_genus')
            self.plot_histogram(metric, df, 'dodge', path / f'{metric}_dodge_hist_no_genus')

    def plot_layered(self, metric: str, df: pd.DataFrame, path: Path):
        g = sns.FacetGrid(df, row='type', height=1.5, aspect=4)
        g.map_dataframe(sns.histplot, x='value', binwidth=self.binwindth, binrange=(0.0, 1.0))
        g.set_axis_labels('distance', 'Count')
        g.set_xlabels(metric)
        for format in self.formats:
            g.savefig(path.with_suffix(f'.{format}'), transparent=True)
        plt.close(g.fig)

    def plot_histogram(self, metric: str, df: pd.DataFrame, multiple: str, path: Path):
        fig, ax = plt.subplots()
        sns.histplot(df, x='value', hue='type', multiple=multiple, binwidth=self.binwindth, binrange=(0.0, 1.0), ax=ax)
        sns.despine()

        ax.margins(x=0.01)
        ax.set_title(ax.get_title(), fontsize=20)
        ax.set_xlabel(ax.get_xlabel(), fontsize=12)
        ax.set_ylabel(ax.get_ylabel(), fontsize=12)
        ax.tick_params(axis='both', labelsize=10)

        for format in self.formats:
            fig.savefig(path.with_suffix(f'.{format}'), transparent=True)
        plt.close(fig)
