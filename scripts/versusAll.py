from __future__ import annotations
from itertools import groupby, product, chain
from pathlib import Path
from sys import argv
from time import perf_counter
from statistics import mean, median, stdev
from typing import NamedTuple, TextIO
from math import inf
from itaxotools.taxi3.align import *
from itaxotools.taxi3.distances import *
from itaxotools.taxi3.pairs import *
from itaxotools.taxi3.sequences import *
from itaxotools.taxi3.partitions import *
from itaxotools.taxi3.handlers import *
from itaxotools.taxi3.statistics import *
from typing import NamedTuple
import numpy as np


class GenericDistance(NamedTuple):
    metric: DistanceMetric
    idx: str
    idy: str
    d: float


class SimpleStatistics(NamedTuple):
    min: float
    max: float
    mean: float


class DistanceStatistics(NamedTuple):
    metric: DistanceMetric
    idx: str
    idy: str
    min: float
    max: float
    mean: float


class SimpleAggregator:
    def __init__(self):
        self.sum = 0.0
        self.min = inf
        self.max = 0.0
        self.n = 0

    def add(self, value: float):
        self.sum += value
        if value < self.min:
            self.min = value
        if value > self.max:
            self.max = value
        self.n += 1

    def calculate(self):
        if not self.n:
            return SimpleStatistics(None, None, None)
        return SimpleStatistics(self.min, self.max, self.sum / self.n)


class DistanceAggregator:
    def __init__(self, metric: DistanceMetric):
        self.metric = metric
        self.aggs = dict()

    def add(self, d: GenericDistance):
        if (d.idx, d.idy) not in self.aggs:
            self.aggs[(d.idx, d.idy)] = SimpleAggregator()
        self.aggs[(d.idx, d.idy)].add(d.d)

    def __iter__(self):
        for (idx, idy), agg in self.aggs.items():
            stats = agg.calculate()
            yield DistanceStatistics(self.metric, idx, idy, stats.min, stats.max, stats.mean)


class SubsetStatisticsHandler(FileHandler[tuple[DistanceStatistics]]):
    def _open(
        self,
        path: Path,
        mode: 'r' | 'w' = 'w',
        missing: str = 'NA',
        formatter: str = '{:f}',
        *args, **kwargs
    ):
        self.missing = missing
        self.formatter = formatter
        super()._open(path, mode, *args, **kwargs)

    def distanceToText(self, d: float | None) -> str:
        if d is None:
            return self.missing
        return self.formatter.format(d)

    def _iter_write(self, *args, **kwargs) -> WriteHandle[tuple[DistanceStatistics]]:
        buffer = None
        with FileHandler.Tabfile(self.path, 'w') as file:
            try:
                bunch = yield
                self._write_headers(file, bunch)
                self._write_stats(file, bunch)
                while True:
                    bunch = yield
                    self._write_stats(file, bunch)
            except GeneratorExit:
                if not buffer:
                    return
                self._write_headers(file, buffer)
                self._write_stats(file, buffer)

    def _write_headers(self, file: TextIO, bunch: tuple[DistanceStatistics]):
        raise NotImplementedError()

    def _write_stats(self, file: TextIO, bunch: tuple[DistanceStatistics]):
        raise NotImplementedError()

    def _iter_read(self, *args, **kwargs) -> None:
        raise NotImplementedError()


class SubsetPairsStatisticsHandler(SubsetStatisticsHandler):
    def _write_headers(self, file: TextIO, bunch: tuple[DistanceStatistics]):
        metrics = (str(stats.metric) for stats in bunch)
        combinations = product(metrics, ['min', 'max', 'mean'])
        headers = (f'{metric} {stat}' for metric, stat in combinations)
        out = ('target', 'query', *headers)
        file.write(out)

    def _write_stats(self, file: TextIO, bunch: tuple[DistanceStatistics]):
        idx = bunch[0].idx
        idy = bunch[0].idy
        stats = ((stats.min, stats.max, stats.mean) for stats in bunch)
        stats = (self.distanceToText(stat) for stat in chain(*stats))
        out = (idx, idy, *stats)
        file.write(out)


class SubsetIdentityStatisticsHandler(SubsetStatisticsHandler):
    def _write_headers(self, file: TextIO, bunch: tuple[DistanceStatistics]):
        metrics = (str(stats.metric) for stats in bunch)
        combinations = product(metrics, ['min', 'max', 'mean'])
        headers = (f'{metric} {stat}' for metric, stat in combinations)
        out = ('target', *headers)
        file.write(out)

    def _write_stats(self, file: TextIO, bunch: tuple[DistanceStatistics]):
        idx = bunch[0].idx
        stats = ((stats.min, stats.max, stats.mean) for stats in bunch)
        stats = (self.distanceToText(stat) for stat in chain(*stats))
        out = (idx, *stats)
        file.write(out)


class SummaryHandler(DistanceHandler.Linear.WithExtras):
    def _open(self, path, mode, species_partition, genera_partition):
        self.species_partition = species_partition
        self.genera_partition = genera_partition
        super()._open(path, mode, tagX=' (query 1)', tagY=' (query 2)')

    def _write_headers(self, file: FileHandler.Tabfile, line: list[Distance]):
        if self.wrote_headers:
            return
        idxHeader = self.idxHeader + self.tagX
        idyHeader = self.idyHeader + self.tagY
        extrasX = [key + self.tagX for key in line[0].x.extras.keys()]
        extrasY = [key + self.tagY for key in line[0].y.extras.keys()]
        metrics = [str(distance.metric) for distance in line]
        infoX = ('genus' + self.tagX, 'species' + self.tagX)
        infoY = ('genus' + self.tagY, 'species' + self.tagY)
        out = (idxHeader, idyHeader, *metrics, *extrasX, *extrasY, *infoX, *infoY, 'comparison_type')
        file.write(out)
        self.wrote_headers = True

    def _write_scores(self, file: FileHandler.Tabfile, line: list[Distance]):
        idx = line[0].x.id
        idy = line[0].y.id
        extrasX = line[0].x.extras.values()
        extrasY = line[0].y.extras.values()
        scores = [self.distanceToText(distance.d) for distance in line]
        genusX = self.genera_partition[idx]
        genusY = self.genera_partition[idy]
        speciesX = self.species_partition[idx]
        speciesY = self.species_partition[idy]
        comparison_type = self._get_comparison_type(genusX, genusY, speciesX, speciesY)
        out = (idx, idy, *scores, *extrasX, *extrasY, genusX, speciesX, genusY, speciesY, comparison_type)
        file.write(out)

    @staticmethod
    def _get_comparison_type(genusX, genusY, speciesX, speciesY) -> str:
        if genusX == genusY:
            if speciesX == speciesY:
                return 'intra-species'
            else:
                return 'inter-species'
        return 'inter-genus'


def multiply(g, n):
    return (x for x in g for i in range(n))


def calculate_statistics_all(data: Sequences, path_stats: Path):
    allStats = StatisticsCalculator()

    for seq in data:
        allStats.add(seq.seq)
        yield seq

    with StatisticsHandler.Single(path_stats, 'w', float_formatter='{:.2f}', percentage_formatter='{:.2f}') as file:
        stats = allStats.calculate()
        file.write(stats)


def calculate_statistics_partition(data: Sequences, partition: Partition, group_name: str, path: Path):
    try:
        calculators = dict()
        for subset in partition.values():
            if subset not in calculators:
                calculators[subset] = StatisticsCalculator(group=subset)

        for seq in data:
            subset = partition[seq.id]
            calculators[subset].add(seq.seq)
            yield seq

    except GeneratorExit:
        pass

    finally:
        with StatisticsHandler.Groups(path, 'w', group_name=group_name, float_formatter='{:.2f}', percentage_formatter='{:.2f}') as file:
            for calc in calculators.values():
                stats = calc.calculate()
                file.write(stats)


def calculate_statistics_species(data: Sequences, partition: Partition, path: Path):
    yield from calculate_statistics_partition(data, partition, 'species', path)


def calculate_statistics_genera(data: Sequences, partition: Partition, path: Path):
    yield from calculate_statistics_partition(data, partition, 'genera', path)


def align_pairs(pairs: SequencePairs):
    aligner = PairwiseAligner.Biopython()
    return aligner.align_pairs(pairs)


def write_pairs(pairs: SequencePairs, path: Path):
    with SequencePairHandler.Formatted(path, 'w') as file:
        for pair in pairs:
            file.write(pair)
            yield pair


def calculate_distances(pairs: SequencePairs, metrics: list[DistanceMetric]):
    for x, y in pairs:
        for metric in metrics:
            yield metric.calculate(x, y)


def write_distances_linear(distances: Distances, path: Path):
    with DistanceHandler.Linear.WithExtras(path, 'w') as file:
        for d in distances:
            file.write(d)
            yield d


def write_distances_matrix(distances: Distances, metric: DistanceMetric, path: Path):
    with DistanceHandler.Matrix(path, 'w') as file:
        for d in distances:
            if d.metric.type == metric.type:
                file.write(d)
            yield d


def write_distances_multimatrix(distances: Distances, metrics: list[DistanceMetric], path: Path):
    for metric in metrics:
        distances = write_distances_matrix(distances, metric, path / f'matrix.{metric.label}.tsv')
    return distances


def aggregate_distances(distances: Distances, partition: Partition, metrics: list[DistanceMetric], group_name: str, path: Path):
    try:
        aggregators = dict()
        for metric in metrics:
            aggregators[str(metric)] = DistanceAggregator(metric)

        for distance in distances:
            subset_x = partition[distance.x.id]
            subset_y = partition[distance.y.id]
            generic = GenericDistance(distance.metric, subset_x, subset_y, distance.d)
            aggregators[str(generic.metric)].add(generic)
            yield generic

    except GeneratorExit:
        pass

    finally:
        with (
            SubsetPairsStatisticsHandler(path / f'{group_name}.pairs.tsv', 'w', formatter='{:.4f}') as pairs_file,
            SubsetIdentityStatisticsHandler(path / f'{group_name}.identity.tsv', 'w', formatter='{:.4f}') as identity_file,
        ):
            aggs = aggregators.values()
            iterators = (iter(agg) for agg in aggs)
            bunches = zip(*iterators)
            for bunch in bunches:
                if bunch[0].idx == bunch[0].idy:
                    identity_file.write(bunch)
                else:
                    pairs_file.write(bunch)


def aggregate_distances_species(distances: Distances, partition: Partition, metrics: list[DistanceMetric], group_name: str, path: Path):
    return aggregate_distances(distances, partition, metrics, group_name, path)


def aggregate_distances_genera(distances: Distances, partition: Partition, metrics: list[DistanceMetric], group_name: str, path: Path):
    return aggregate_distances(distances, partition, metrics, group_name, path)


def write_summary(distances, species_partition, genera_partition, path: Path):
    with SummaryHandler(path, 'w', species_partition, genera_partition) as file:
        for distance in distances:
            if 'organism' in distance.x.extras:
                del distance.x.extras['organism']
            if 'organism' in distance.y.extras:
                del distance.y.extras['organism']
            file.write(distance)
            yield distance


def progress(distances, total):
    for index, distance in enumerate(distances, 1):
        print(f"\rCalculating... {index}/{total} = {100*index/total:.2f}%", end="")
        yield distance
    print('\nFinalizing...')


def main():
    path_data = Path(argv[1])
    path_out = Path(argv[2])
    path_out.mkdir(exist_ok=True)
    path_stats = path_out / 'stats'
    path_stats.mkdir(exist_ok=True)
    path_align = path_out / 'align'
    path_align.mkdir(exist_ok=True)
    path_distances = path_out / 'distances'
    path_distances.mkdir(exist_ok=True)
    path_subs = path_out / 'subsets'
    path_subs.mkdir(exist_ok=True)

    metrics = [
        DistanceMetric.Uncorrected(),
        DistanceMetric.UncorrectedWithGaps(),
        DistanceMetric.JukesCantor(),
        DistanceMetric.Kimura2P(),
    ]

    ts = perf_counter()

    species_partition = Partition.fromPath(path_data, PartitionHandler.Tabfile, idHeader='seqid', subHeader='organism')
    genera_partition = Partition.fromPath(path_data, PartitionHandler.Tabfile, idHeader='seqid', subHeader='organism', filter=PartitionHandler.subset_first_word)

    data = Sequences.fromPath(path_data, SequenceHandler.Tabfile, idHeader='seqid', seqHeader='sequence')
    data = data.normalize()

    data_left = data
    data_left = calculate_statistics_all(data_left, path_stats / 'all.tsv')
    data_left = calculate_statistics_species(data_left, species_partition, path_stats / 'species.tsv')
    data_left = calculate_statistics_genera(data_left, genera_partition, path_stats / 'genera.tsv')

    pairs = SequencePairs.fromProduct(data_left, data)
    pairs = align_pairs(pairs)
    pairs = write_pairs(pairs, path_align / 'alignments.txt')

    distances = calculate_distances(pairs, metrics)
    distances = write_distances_linear(distances, path_distances / 'linear.tsv')
    distances = write_distances_multimatrix(distances, metrics, path_distances)

    distances = multiply(distances, 3)
    distances_species = aggregate_distances_species(distances, species_partition, metrics, 'species', path_subs)
    distances_genera = aggregate_distances_genera(distances, genera_partition, metrics, 'genera', path_subs)
    distances = (distances for distances, _, _ in zip(distances, distances_species, distances_genera))

    distances = progress(distances, len(metrics) * len(data) ** 2)

    distances = write_summary(distances, species_partition, genera_partition, path_out / 'summary.tsv')

    for x in distances:
        pass

    tf = perf_counter()
    print(f'Time taken: {tf-ts:.4f}s')

    return


if __name__ == '__main__':
    main()
