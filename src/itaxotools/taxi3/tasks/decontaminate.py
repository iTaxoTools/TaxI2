from __future__ import annotations

from itertools import groupby, product, chain
from pathlib import Path
from time import perf_counter
from statistics import mean, median, stdev
from typing import NamedTuple, TextIO, Callable
from math import inf

from itaxotools.common.utility import AttrDict

from ..align import PairwiseAligner
from ..distances import Distances, DistanceHandler, DistanceMetric
from ..pairs import SequencePair, SequencePairs, SequencePairHandler
from ..sequences import Sequence, Sequences, SequenceHandler
from ..partitions import Partition, PartitionHandler
from ..statistics import StatisticsCalculator, StatisticsHandler
from ..handlers import FileHandler, ReadHandle, WriteHandle
from ..plot import HistogramPlotter


def multiply(iterator: iter, n: int):
    return (item for item in iterator for i in range(n))


def split(source: iter, *funcs: list[Callable]):
    source = multiply(source, len(funcs))
    return [
        map(func, source)
        for func in funcs
    ]


def console_report(caption, index, total):
    if caption == 'Finalizing...':
        print(f"\rCalculating... {total}/{total} = {100:.2f}%", end="")
        print('\nFinalizing...')
    else:
        print(f"\rCalculating... {index}/{total} = {100*index/total:.2f}%", end="")


class Verdict(NamedTuple):
    sequence: Sequence
    contaminant: bool


class SummaryLine(NamedTuple):
    query_id: str
    outgroup_id: str
    outgroup_distance: float
    contaminant: bool


class Results(NamedTuple):
    output_directory: Path
    seconds_taken: float


class SummaryHandle(FileHandler[SummaryLine]):
    def _open(
        self,
        path: Path,
        mode: 'r' | 'w' = 'w',
        formatter: str = '{:f}',
        *args, **kwargs
    ):
        self.formatter = formatter
        super()._open(path, mode, *args, **kwargs)

    def _iter_read(self, *args, **kwargs) -> ReadHandle[SummaryLine]:
        raise NotImplementedError()

    def format_line(self, line: SummaryLine):
        return (
            line.query_id,
            line.outgroup_id,
            self.formatter.format(line.outgroup_distance),
            'Yes' if line.contaminant else 'No',
        )

    def _iter_write(self, *args, **kwargs) -> WriteHandle[SummaryLine]:
        try:
            with FileHandler.Tabfile(self.path, 'w') as file:
                while True:
                    line = yield
                    file.write(self.format_line(line))
        except GeneratorExit:
            return


class Decontaminate:

    def __init__(self):

        self.work_dir: Path = None
        self.paths = AttrDict()

        self.progress_handler: Callable = console_report
        self.progress_interval: float = 0.015

        self.input: Sequences = None
        self.outgroup: Sequences = None

        self.params = AttrDict()

        self.params.thresholds = AttrDict()
        self.params.thresholds.similarity: float = 0.07

        self.params.pairs = AttrDict()
        self.params.pairs.align: bool = True
        self.params.pairs.write: bool = True
        self.params.pairs.scores: Scores = None

        self.params.distances = AttrDict()
        self.params.distances.metric: DistanceMetric = None
        self.params.distances.write_linear: bool = True
        self.params.distances.write_matricial: bool = True

        self.params.format = AttrDict()
        self.params.format.float: str = '{:.4f}'
        self.params.format.missing: str = 'NA'
        self.params.format.percentage_multiply: bool = False

    def check_metric(self):
        self.params.distances.metric = self.params.distances.metric or DistanceMetric.Uncorrected()

    def generate_paths(self):
        assert self.work_dir
        self.create_parents(self.work_dir)
        metric = str(self.params.distances.metric)
        extension = '.tsv'

        self.paths.summary = self.work_dir / 'summary.tsv'
        self.paths.decontaminated = self.work_dir / f'decontaminated{extension}'
        self.paths.contaminants = self.work_dir / f'contaminants{extension}'
        self.paths.aligned_pairs = self.work_dir / 'aligned_pairs.txt'
        self.paths.distances_linear = self.work_dir / 'distances' / f'{metric}.linear.tsv'
        self.paths.distances_matrix = self.work_dir / 'distances' / f'{metric}.matricial.tsv'

    def create_parents(self, path: Path):
        if path.suffix:
            path = path.parent
        path.mkdir(parents=True, exist_ok=True)

    def align_pairs(self, pairs: iter[SequencePair]) -> iter[SequencePair]:
        if not self.params.pairs.align:
            yield from pairs
            return

        aligner = PairwiseAligner.Biopython(self.params.pairs.scores)
        yield from aligner.align_pairs(pairs)

    def write_pairs(self, pairs: iter[SequencePair]) -> iter[SequencePair]:
        if not self.params.pairs.write:
            yield from pairs
            return

        self.create_parents(self.paths.aligned_pairs)
        with SequencePairHandler.Formatted(self.paths.aligned_pairs, 'w') as file:
            for pair in pairs:
                file.write(pair)
                yield pair

    def calculate_distances(self, pairs: iter[SequencePair]) -> iter[Distance]:
        metric = self.params.distances.metric
        for x, y in pairs:
            yield metric.calculate(x, y)

    def adjust_distances(self, distances: iter[Distance]) -> iter[Distance]:
        if not self.params.format.percentage_multiply:
            yield from distances
            return

        for distance in distances:
            distance = distance._replace(d = distance.d * 100)
            yield distance

    def write_distances_linear(self, distances: iter[Distance], path: Path) -> iter[Distance]:
        if not self.params.distances.write_linear:
            yield from distances
            return

        self.create_parents(path)
        with DistanceHandler.Linear.WithExtras(
            path, 'w',
            missing = self.params.format.missing,
            formatter = self.params.format.float,
        ) as file:
            for distance in distances:
                file.write(distance)
                yield distance

    def write_outgroup_distances_linear(self, distances: iter[Distance]) -> iter[Distance]:
        return self.write_distances_linear(distances, self.paths.distances_linear)

    def write_distances_matrix(self, distances: iter[Distance], path: Path) -> iter[Distance]:
        if not self.params.distances.write_matricial:
            yield from distances
            return

        self.create_parents(path)
        with DistanceHandler.Matrix(
            path, 'w',
            missing = self.params.format.missing,
            formatter = self.params.format.float,
        ) as file:
            for distance in distances:
                file.write(distance)
                yield distance

    def write_outgroup_distances_matrix(self, distances: iter[Distance]) -> iter[Distance]:
        return self.write_distances_matrix(distances, self.paths.distances_matrix)

    def group_distances_left(self, distances: iter[Distance]) -> iter[iter[Distance]]:
        for _, group in groupby(distances, lambda distance: distance.x.id):
            yield group

    def get_minimum_distances(self, groups: iter[iter[Distance]]) -> iter[Distance]:
        for distances in groups:
            distances = (distance for distance in distances if distance.d is not None)
            minimum_distance = min(distances, key = lambda distance: distance.d)
            yield minimum_distance

    def _find_contaminants(
        self,
        sequences: iter[Sequence],
        out_minimums: iter[Distance],
    ) -> iter[tuple[Verdict, SummaryLine]]:

        threshold = self.params.thresholds.similarity
        all = zip(sequences, out_minimums)
        for sequence, outgroup_minimum in all:
            is_contaminant = bool(outgroup_minimum.d <= threshold)
            verdict = Verdict(sequence, is_contaminant)
            line = SummaryLine(
                query_id = sequence.id,
                outgroup_id = outgroup_minimum.y.id,
                outgroup_distance = outgroup_minimum.d,
                contaminant = is_contaminant,
            )
            yield (verdict, line)

    def find_contaminants(
        self,
        sequences: iter[Sequence],
        out_minimums: iter[Distance],
    ) -> tuple[iter[Verdict], iter[SummaryLine]]:

        data = self._find_contaminants(sequences, out_minimums)
        return split(data, lambda x: x[0], lambda x: x[1])

    def write_file_decontaminated(self, verdicts: iter[Verdict]) -> iter[Verdict]:
        with SequenceHandler.Tabfile(self.paths.decontaminated, 'w') as file:
            for verdict in verdicts:
                if not verdict.contaminant:
                    file.write(verdict.sequence)
                yield verdict

    def write_file_contaminants(self, verdicts: iter[Verdict]) -> iter[Verdict]:
        with SequenceHandler.Tabfile(self.paths.contaminants, 'w') as file:
            for verdict in verdicts:
                if verdict.contaminant:
                    file.write(verdict.sequence)
                yield verdict

    def write_summary(self, lines: iter[SummaryLine]) -> iter[SummaryLine]:
        with SummaryHandle(
            self.paths.summary, 'w',
            formatter = self.params.format.float,
        ) as file:
            for line in lines:
                if line is not None:
                    file.write(line)
                yield line

    def report_progress(self, verdicts: iter[Verdict]):
        total = len(self.input)
        last_time = perf_counter()
        for index, verdict in enumerate(verdicts, 1):
            new_time = perf_counter()
            if new_time - last_time >= self.progress_interval:
                self.progress_handler('verdict.x.id', index, total)
                last_time = new_time
            yield verdict
        self.progress_handler('Finalizing...', total, total)

    def start(self) -> None:
        ts = perf_counter()

        self.check_metric()
        self.generate_paths()

        data = self.input.normalize()

        outgroup = self.outgroup.normalize()
        out_pairs = SequencePairs.fromProduct(data, outgroup)
        out_pairs = self.align_pairs(out_pairs)
        out_pairs = self.write_pairs(out_pairs)
        out_distances = self.calculate_distances(out_pairs)
        out_distances = self.adjust_distances(out_distances)
        out_distances = self.write_outgroup_distances_linear(out_distances)
        out_distances = self.write_outgroup_distances_matrix(out_distances)
        out_groups = self.group_distances_left(out_distances)
        out_minimums = self.get_minimum_distances(out_groups)

        verdicts, lines = self.find_contaminants(data, out_minimums)
        verdicts = self.write_file_decontaminated(verdicts)
        verdicts = self.write_file_contaminants(verdicts)
        verdicts = self.report_progress(verdicts)
        lines = self.write_summary(lines)

        for _ in zip(verdicts, lines):
            pass

        tf = perf_counter()

        return Results(self.work_dir, tf - ts)
