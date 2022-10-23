from __future__ import annotations

from typing import Callable, Iterable, NamedTuple
from enum import Enum, auto

from .sequences import Sequence, Sequences
from .types import Type, Container


class Distance(NamedTuple):
    metric: DistanceMetric
    idx: str
    idy: str
    d: float | None


class Distances(Container[Distance]):
    @classmethod
    def fromFile(cls, file: DistanceFile, *args, **kwargs) -> Distances:
        return cls(file.read, *args, **kwargs)


class DistanceFile(Type):
    """Handlers for distance files"""

    def __init__(self, path: Path):
        self.path = path

    def read(self, *args, **kwargs) -> iter[Distance]:
        raise NotImplementedError()

    def write(self, distances: iter[Distance], *args, **kwargs) -> None:
        raise NotImplementedError()


class Buffer(DistanceFile):
    def __init__(self):
        self.distances = []

    def read(self) -> iter[Distance]:
        return iter(self.distances)

    def write(self, distances: iter[Distance]) -> None:
        for distance in distances:
            self.distances.append(distance)


class Linear(DistanceFile):
    def read(self) -> iter[Distance]:
        with open(self.path, 'r') as f:
            data = f.readline()
            id1Header, id2Header, label = data.strip().split('\t')
            metric = DistanceMetric.fromLabel(label)
            for line in f:
                id1, id2, labelData = line.split('\t')
                yield Distance(metric, id1, id2, float(labelData))


class Matrix(DistanceFile):
    def read(self, metric: DistanceMetric = None) -> iter[Distance]:
        metric = metric or DistanceMetric.Unknown()
        raise NotImplementedError()


class DistanceMetric(Type):
    """Metrics for calculating distances"""
    label: str

    def __str__(self):
        return self.label

    def _calculate(self, x: str, y: str) -> float:
        # Legacy implementation found in these modules:
        # library.task.CalculateDistances
        # library.datatypes.Metric
        # library.rust_backend
        # library.calculate_distances
        # library.alfpy_distance
        # alfpy.bbc
        # alfpy.ncd
        # Unlike legacy code, we need to calculate distances between
        # one pair of sequences at a time. Bulk processing comes later.
        # Must dig in the code and find the simplest solution.
        raise NotImplementedError()

    def calculate(self, x: Sequence, y: Sequence) -> Distance:
        return Distance(self, x.id, y.id, self._calculate(x.seq, y.seq))

    @classmethod
    def fromLabel(cls, label: str):
        for child in cls:
            if label == child.label:
                return child()


class Unknown(DistanceMetric):
    label = '?'


class Uncorrected(DistanceMetric):
    label = 'p-distance'


class UncorrectedWithGaps(DistanceMetric):
    label = 'p-distance with gaps'


class JukesCantor(DistanceMetric):
    label = 'jc'


class Kimura2P(DistanceMetric):
    label = 'k2p'


class NCD(DistanceMetric):
    label = 'ncd({})'

    def __init__(self, arg):
        self.arg = arg

    def __str__(self):
        return self.label.format(self.arg)

    def __eq__(self, other):
        return super().__eq__(other) and self.arg == other.arg


class BBC(DistanceMetric):
    label = 'bbc({})'

    def __init__(self, arg):
        self.arg = arg

    def __str__(self):
        return self.label.format(self.arg)

    def __eq__(self, other):
        return super().__eq__(other) and self.arg == other.arg
