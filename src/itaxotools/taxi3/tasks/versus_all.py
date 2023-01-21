from __future__ import annotations

from itertools import groupby, product, chain
from pathlib import Path
from time import perf_counter
from statistics import mean, median, stdev
from typing import NamedTuple, TextIO
from math import inf

from itaxotools.common.utility import AttrDict

from ..align import PairwiseAligner
from ..distances import Distances, DistanceHandler, DistanceMetric
from ..pairs import SequencePairs, SequencePairHandler
from ..sequences import Sequences, SequenceHandler
from ..partitions import Partition, PartitionHandler
from ..statistics import StatisticsCalculator, StatisticsHandler
from ..handlers import FileHandler


def multiply(iterator: iter, n: int):
    return (item for item in iterator for i in range(n))


def console_report(caption, index, total):
    if caption == 'Finalizing...':
        print('\nFinalizing...')
    else:
        print(f"\rCalculating... {index}/{total} = {100*index/total:.2f}%", end="")


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


class SubsetPair(NamedTuple):
    x: str
    y: str


class SubsetDistance(NamedTuple):
    distance: Distance
    genera: SubsetPair | None
    species: SubsetPair | None

    def get_comparison_type(self) -> str:
        same_genera = bool(self.genera.x == self.genera.y) if self.genera else None
        same_species = bool(self.species.x == self.species.y) if self.species else None
        return {
            (None, None): '-',
            (None, True): 'intra-species',
            (None, False): 'inter-species',
            (False, None): 'inter-genus',
            (False, True): 'inter-genus',
            (False, False): 'inter-genus',
            (True, None): 'intra-genus',
            (True, True): 'intra-species',
            (True, False): 'inter-species',
        }[(same_genera, same_species)]



class SummaryHandler(DistanceHandler.Linear.WithExtras):
    def _open(self, path, mode):
        super()._open(path, mode, tagX=' (query 1)', tagY=' (query 2)')

    def _assemble_line(self) -> Generator[None, SubsetDistance, list[SubsetDistance]]:
        buffer = self.buffer
        try:
            while True:
                distance = yield
                buffer.append(distance)
                if any((
                    buffer[0].distance.x.id != buffer[-1].distance.x.id,
                    buffer[0].distance.y.id != buffer[-1].distance.y.id,
                )):
                    self.buffer = buffer[-1:]
                    return buffer[:-1]
        except GeneratorExit:
            return

    def _write_headers(self, file: FileHandler.Tabfile, line: list[SubsetDistance]):
        if self.wrote_headers:
            return
        idxHeader = self.idxHeader + self.tagX
        idyHeader = self.idyHeader + self.tagY
        extrasX = [key + self.tagX for key in line[0].distance.x.extras.keys()]
        extrasY = [key + self.tagY for key in line[0].distance.y.extras.keys()]
        metrics = [str(subset_distance.distance.metric) for subset_distance in line]
        infoX = ('genus' + self.tagX, 'species' + self.tagX)
        infoY = ('genus' + self.tagY, 'species' + self.tagY)
        out = (idxHeader, idyHeader, *metrics, *extrasX, *extrasY, *infoX, *infoY, 'comparison_type')
        file.write(out)
        self.wrote_headers = True

    def _write_scores(self, file: FileHandler.Tabfile, line: list[SubsetDistance]):
        first = line[0]
        idx = first.distance.x.id
        idy = first.distance.y.id
        extrasX = first.distance.x.extras.values()
        extrasY = first.distance.y.extras.values()
        scores = [self.distanceToText(subset_distance.distance.d) for subset_distance in line]
        genusX = first.genera.x if first.genera else '-'
        genusY = first.genera.y if first.genera else '-'
        speciesX = first.species.x if first.species else '-'
        speciesY = first.species.y if first.species else '-'
        comparison_type = first.get_comparison_type()
        out = (idx, idy, *scores, *extrasX, *extrasY, genusX or '-', speciesX or '-', genusY or '-', speciesY or '-', comparison_type)
        file.write(out)


class Results(NamedTuple):
    output_directory: Path
    seconds_taken: float

#     stats_all: Path | None
#     stats_species: Path | None
#     stats_genus: Path | None
#
#     aligned_pairs: Path | None
#
#     summary: Path
#
#     ALOTOFFILES: ...
#
#     matrices: dict[str, Path]
#
#
#     number_of_files_created: int


class VersusAll:

    def __init__(self):

        self.work_dir: Path = None
        self.paths = AttrDict()

        self.progress_handler: Callable = console_report

        self.input = AttrDict()
        self.input.sequences: Sequences = None
        self.input.species: Partition = None
        self.input.genera: Partition = None

        self.params = AttrDict()

        self.params.pairs = AttrDict()
        self.params.pairs.align: bool = True
        self.params.pairs.write: bool = True
        self.params.pairs.scores: Scores = None

        self.params.distances = AttrDict()
        self.params.distances.metrics: list[DistanceMetric] = None
        self.params.distances.write_linear: bool = True
        self.params.distances.write_matricial: bool = True

        self.params.format = AttrDict()
        self.params.format.float: str = '{:.4f}'
        self.params.format.percentage: str = '{:.2f}'
        self.params.format.missing: str = 'NA'
        self.params.format.percentage_multiply: bool = False

        self.params.stats = AttrDict()
        self.params.stats.all: bool = True
        self.params.stats.species: bool = True
        self.params.stats.genera: bool = True

    def generate_paths(self):
        assert self.work_dir

        self.paths.summary = self.work_dir / 'summary.tsv'
        self.paths.stats_all = self.work_dir / 'stats' / 'all.tsv'
        self.paths.stats_species = self.work_dir / 'stats' / 'species.tsv'
        self.paths.stats_genera = self.work_dir / 'stats' / 'genera.tsv'
        self.paths.aligned_pairs = self.work_dir / 'align' / 'aligned_pairs.txt'
        self.paths.distances_linear = self.work_dir / 'distances' / 'linear.tsv'
        self.paths.distances_matricial = self.work_dir / 'distances' / 'matricial'
        self.paths.subsets = self.work_dir / 'subsets'

        for path in [
            self.paths.summary,
        ]:
            self.create_parents(path)

    def create_parents(self, path: Path):
        if path.suffix:
            path = path.parent
        path.mkdir(parents=True, exist_ok=True)

    def check_metrics(self):
        self.params.distances.metrics = self.params.distances.metrics or [
                DistanceMetric.Uncorrected(),
                DistanceMetric.UncorrectedWithGaps(),
                DistanceMetric.JukesCantor(),
                DistanceMetric.Kimura2P(),
            ]

    def calculate_statistics_all(self, sequences: Sequences):
        if not self.params.stats.all:
            yield from sequences
            return

        allStats = StatisticsCalculator()

        for sequence in sequences:
            allStats.add(sequence.seq)
            yield sequence

        self.create_parents(self.paths.stats_all)
        with StatisticsHandler.Single(self.paths.stats_all, 'w', float_formatter='{:.2f}', percentage_formatter='{:.2f}') as file:
            stats = allStats.calculate()
            file.write(stats)

    def calculate_statistics_species(self, sequences: Sequences):
        if not self.input.species:
            return sequences
        if not self.params.stats.species:
            return sequences

        return self._calculate_statistics_partition(sequences, self.input.species, 'species', self.paths.stats_species)

    def calculate_statistics_genera(self, sequences: Sequences):
        if not self.input.genera:
            return sequences
        if not self.params.stats.genera:
            return sequences

        return self._calculate_statistics_partition(sequences, self.input.genera, 'genera', self.paths.stats_genera)

    def _calculate_statistics_partition(self, sequences: Sequences, partition: Partition, group_name: str, path: Path):
        try:
            calculators = dict()
            for subset in partition.values():
                if subset not in calculators:
                    calculators[subset] = StatisticsCalculator(group=subset)

            for sequence in sequences:
                subset = partition[sequence.id]
                calculators[subset].add(sequence.seq)
                yield sequence

        except GeneratorExit:
            pass

        finally:
            self.create_parents(path)
            with StatisticsHandler.Groups(path, 'w', group_name=group_name, float_formatter='{:.2f}', percentage_formatter='{:.2f}') as file:
                for calc in calculators.values():
                    stats = calc.calculate()
                    file.write(stats)

    def align_pairs(self, pairs: SequencePairs):
        if not self.params.pairs.align:
            yield from pairs
            return

        aligner = PairwiseAligner.Biopython()
        yield from aligner.align_pairs(pairs)

    def write_pairs(self, pairs: SequencePairs):
        if not self.params.pairs.write:
            yield from pairs
            return

        self.create_parents(self.paths.aligned_pairs)
        with SequencePairHandler.Formatted(self.paths.aligned_pairs, 'w') as file:
            for pair in pairs:
                file.write(pair)
                yield pair

    def calculate_distances(self, pairs: SequencePairs):
        for x, y in pairs:
            for metric in self.params.distances.metrics:
                yield metric.calculate(x, y)

    def write_distances_linear(self, distances: Distances):
        if not self.params.distances.write_linear:
            yield from distances
            return

        self.create_parents(self.paths.distances_linear)
        with DistanceHandler.Linear.WithExtras(self.paths.distances_linear, 'w') as file:
            for distance in distances:
                file.write(distance)
                yield distance

    def write_distances_multimatrix(self, distances: Distances):
        if not self.params.distances.write_matricial:
            return distances

        self.create_parents(self.paths.distances_matricial)
        for metric in self.params.distances.metrics:
            distances = self._write_distances_matrix(distances, metric, self.paths.distances_matricial / f'{metric.label}.tsv')
        return distances

    def _write_distances_matrix(self, distances: Distances, metric: DistanceMetric, path: Path):
        with DistanceHandler.Matrix(path, 'w') as file:
            for distance in distances:
                if distance.metric.type == metric.type:
                    file.write(distance)
                yield distance

    def aggregate_distances_species(self, distances: Distances) -> iter[SubsetPair | None]:
        if not self.input.species:
            return (None for _ in distances)
        return self._aggregate_distances(distances, self.input.species, 'species', self.paths.subsets)

    def aggregate_distances_genera(self, distances: Distances) -> iter[SubsetPair | None]:
        if not self.input.genera:
            return (None for _ in distances)
        return self._aggregate_distances(distances, self.input.genera, 'genera', self.paths.subsets)

    def _aggregate_distances(self, distances: Distances, partition: Partition, group_name: str, path: Path) -> iter[SubsetPair]:
        try:
            aggregators = dict()
            for metric in self.params.distances.metrics:
                aggregators[str(metric)] = DistanceAggregator(metric)

            for distance in distances:
                subset_x = partition[distance.x.id]
                subset_y = partition[distance.y.id]
                generic = GenericDistance(distance.metric, subset_x, subset_y, distance.d)
                aggregators[str(generic.metric)].add(generic)
                yield SubsetPair(subset_x, subset_y)

        except GeneratorExit:
            pass

        finally:
            self.create_parents(path)
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

    def write_summary(self, distances: iter[SubsetDistance]):
        with SummaryHandler(self.paths.summary, 'w') as file:
            for distance in distances:
                # if 'organism' in distance.x.extras:
                #     del distance.x.extras['organism']
                # if 'organism' in distance.y.extras:
                #     del distance.y.extras['organism']
                file.write(distance)
                yield distance

    def report_progress(self, distances: iter[SubsetDistance]):
        total = len(self.params.distances.metrics) * len(self.input.sequences) ** 2
        from sys import stderr
        for index, distance in enumerate(distances, 1):
            self.progress_handler('distance.x.id', index, total)
            print(index, str(distance.distance.metric), distance.distance.x.id, distance.distance.y.id, file=stderr)
            yield distance
        self.progress_handler('Finalizing...', 0, 0)

    def start(self) -> None:
        ts = perf_counter()

        self.generate_paths()
        self.check_metrics()

        sequences = self.input.sequences.normalize()

        sequences_left = sequences
        sequences_left = self.calculate_statistics_all(sequences_left)
        sequences_left = self.calculate_statistics_species(sequences_left)
        sequences_left = self.calculate_statistics_genera(sequences_left)

        pairs = SequencePairs.fromProduct(sequences_left, sequences)
        pairs = self.align_pairs(pairs)
        pairs = self.write_pairs(pairs)

        distances = self.calculate_distances(pairs)
        distances = self.write_distances_linear(distances)
        distances = self.write_distances_multimatrix(distances)

        distances = multiply(distances, 3)
        genera_pair = self.aggregate_distances_genera(distances)
        species_pair = self.aggregate_distances_species(distances)
        subset_distances = (SubsetDistance(d, g, s) for d, g, s in zip(distances, genera_pair, species_pair))

        subset_distances = self.write_summary(subset_distances)

        subset_distances = self.report_progress(subset_distances)

        for _ in subset_distances:
            pass

        tf = perf_counter()

        return Results(self.work_dir, tf - ts)
