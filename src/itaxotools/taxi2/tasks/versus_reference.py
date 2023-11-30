from __future__ import annotations

from itertools import groupby
from pathlib import Path
from time import perf_counter
from typing import Callable, NamedTuple

from itaxotools.common.utility import AttrDict

from ..align import PairwiseAligner, Scores
from ..distances import Distance, DistanceHandler, DistanceMetric, Distances
from ..pairs import SequencePairHandler, SequencePairs
from ..sequences import Sequences


def multiply(iterator: iter, n: int):
    return (item for item in iterator for i in range(n))


def console_report(caption, index, total):
    if caption == "Finalizing...":
        print(f"\rCalculating... {total}/{total} = {100:.2f}%", end="")
        print("\nFinalizing...")
    else:
        print(f"\rCalculating... {index}/{total} = {100*index/total:.2f}%", end="")


class Results(NamedTuple):
    output_directory: Path
    seconds_taken: float


class VersusReference:
    def __init__(self):
        self.work_dir: Path = None
        self.paths = AttrDict()

        self.progress_handler: Callable = console_report
        self.progress_interval: float = 0.015

        self.input = AttrDict()
        self.input.data: Sequences = None
        self.input.reference: Sequences = None

        self.params = AttrDict()

        self.params.pairs = AttrDict()
        self.params.pairs.align: bool = True
        self.params.pairs.write: bool = True
        self.params.pairs.scores: Scores = None

        self.params.distances = AttrDict()
        self.params.distances.metric: DistanceMetric = None
        self.params.distances.extra_metrics: list[DistanceMetric] = None
        self.params.distances.write_linear: bool = True
        self.params.distances.write_matricial: bool = True

        self.params.format = AttrDict()
        self.params.format.float: str = "{:.4f}"
        self.params.format.percentage: str = "{:.2f}"
        self.params.format.missing: str = "NA"
        self.params.format.percentage_multiply: bool = False

    def generate_paths(self):
        assert self.work_dir
        self.create_parents(self.work_dir)
        metric = str(self.params.distances.metric)

        self.paths.closest = self.work_dir / "closest.tsv"
        self.paths.aligned_pairs = self.work_dir / "aligned_pairs.txt"
        self.paths.distances_linear = (
            self.work_dir / "distances" / f"{metric}.linear.tsv"
        )
        self.paths.distances_matricial = (
            self.work_dir / "distances" / f"{metric}.matricial.tsv"
        )

    def create_parents(self, path: Path):
        if path.suffix:
            path = path.parent
        path.mkdir(parents=True, exist_ok=True)

    def check_metrics(self):
        self.params.distances.metric = (
            self.params.distances.metric or DistanceMetric.Uncorrected()
        )
        self.params.distances.extra_metrics = self.params.distances.extra_metrics or [
            DistanceMetric.UncorrectedWithGaps(),
            DistanceMetric.JukesCantor(),
            DistanceMetric.Kimura2P(),
        ]
        if self.params.distances.metric in self.params.distances.extra_metrics:
            self.params.distances.extra_metrics.remove(self.params.distances.metric)

    def normalize_sequences(self, sequences: Sequences) -> Sequences:
        if not self.params.pairs.align:
            return sequences
        return sequences.normalize()

    def align_pairs(self, pairs: SequencePairs):
        if not self.params.pairs.align:
            yield from pairs
            return

        aligner = PairwiseAligner.Biopython(self.params.pairs.scores)
        yield from aligner.align_pairs(pairs)

    def write_pairs(self, pairs: SequencePairs):
        if not self.params.pairs.write:
            yield from pairs
            return

        self.create_parents(self.paths.aligned_pairs)
        with SequencePairHandler.Formatted(self.paths.aligned_pairs, "w") as file:
            for pair in pairs:
                file.write(pair)
                yield pair

    def calculate_distances(self, pairs: SequencePairs):
        metric = self.params.distances.metric
        for x, y in pairs:
            yield metric.calculate(x, y)

    def calculate_extra_distances(self, distances: Distances):
        metrics = self.params.distances.extra_metrics
        for distance in distances:
            yield distance
            for metric in metrics:
                yield metric.calculate(distance.x, distance.y)

    def adjust_distances(self, distances: Distances):
        if not self.params.format.percentage_multiply:
            yield from distances
            return

        for distance in distances:
            if distance.d is not None:
                distance = distance._replace(d=distance.d * 100)
            yield distance

    def adjust_extra_distances(self, distances: Distances):
        if not self.params.format.percentage_multiply:
            yield from distances
            return

        for distance in distances:
            if distance.metric != self.params.distances.metric:
                if distance.d is not None:
                    distance = distance._replace(d=distance.d * 100)
            yield distance

    def write_distances_linear(self, distances: Distances):
        if not self.params.distances.write_linear:
            yield from distances
            return

        self.create_parents(self.paths.distances_linear)
        with DistanceHandler.Linear.WithExtras(
            self.paths.distances_linear,
            "w",
            missing=self.params.format.missing,
            formatter=self.params.format.float,
        ) as file:
            for distance in distances:
                file.write(distance)
                yield distance

    def write_distances_matrix(self, distances: Distances):
        if not self.params.distances.write_matricial:
            yield from distances
            return

        self.create_parents(self.paths.distances_matricial)
        with DistanceHandler.Matrix(
            self.paths.distances_matricial,
            "w",
            missing=self.params.format.missing,
            formatter=self.params.format.float,
        ) as file:
            for distance in distances:
                file.write(distance)
                yield distance

    def get_minimum_distances(self, distances: Distances):
        for _, group in groupby(distances, lambda distance: distance.x.id):
            group = (distance for distance in group if distance.d is not None)
            minimum_distance = min(group, key=lambda distance: distance.d)
            yield minimum_distance

    def write_closest_distances(self, distances: iter[Distance]):
        self.create_parents(self.paths.closest)
        with DistanceHandler.Linear.WithExtras(
            self.paths.closest,
            "w",
            missing=self.params.format.missing,
            formatter=self.params.format.float,
        ) as file:
            for distance in distances:
                file.write(distance)
                yield distance

    def report_progress(self, distances: iter[Distance]):
        total = len(self.input.data) * len(self.input.reference)
        last_time = perf_counter()
        for index, distance in enumerate(distances, 1):
            new_time = perf_counter()
            if new_time - last_time >= self.progress_interval:
                self.progress_handler("distance.x.id", index, total)
                last_time = new_time
            yield distance
        self.progress_handler("Finalizing...", total, total)

    def start(self) -> None:
        ts = perf_counter()

        self.check_metrics()
        self.generate_paths()

        data = self.input.data
        data = self.normalize_sequences(data)

        reference = self.input.reference
        reference = self.normalize_sequences(reference)

        pairs = SequencePairs.fromProduct(data, reference)
        pairs = self.align_pairs(pairs)
        pairs = self.write_pairs(pairs)

        distances = self.calculate_distances(pairs)
        distances = self.report_progress(distances)

        distances = self.adjust_distances(distances)
        distances = self.write_distances_linear(distances)
        distances = self.write_distances_matrix(distances)

        minimum_distances = self.get_minimum_distances(distances)

        all_distances = self.calculate_extra_distances(minimum_distances)
        all_distances = self.adjust_extra_distances(all_distances)
        all_distances = self.write_closest_distances(all_distances)

        for _ in all_distances:
            pass

        tf = perf_counter()

        return Results(self.work_dir, tf - ts)
