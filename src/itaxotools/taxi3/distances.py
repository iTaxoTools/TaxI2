from __future__ import annotations

from typing import Callable, NamedTuple
from enum import Enum, auto

from .sequences import Sequence, Sequences
from .types import Type


class Distance(NamedTuple):
    metric: DistanceMetric
    idx: str
    idy: str
    d: float


class Distances:
    """Distance containers that can be iterated multiple times"""
    def __init__(self, source: Callable[None, iter[Distance]]):
        self.source = source

    def __iter__(self):
        return self.source()

    @classmethod
    def fromFile(cls, file: DistanceFile):
        return cls(file.read)


class DistanceFile(Type):
    """Handlers for distance files"""
    def __init__(self, path: Path):
        self.path = path

    def read(self, **kwargs) -> iter[Distance]:
        raise NotImplementedError()

    def write(self, distances: iter[Distance]) -> None:
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
    pass


class Matrix(DistanceFile):
    pass


class DistanceMetric(Type):
    """Metrics for calculating distances"""
    label: str

    def __str__(self):
        return self.label

    def _calculate(self, x: str, y: str) -> float:
        raise NotImplementedError()

    def calculate(self, x: Sequence, y: Sequence) -> Distance:
        return Distance(self, x.id, y.id, self._calculate(x.seq, y.seq))


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

    def __str__(self):
        return self.label.format(self.arg)

    def __init__(self, arg):
        self.arg = arg


class BBC(DistanceMetric):
    label = 'bbc({})'

    def __str__(self):
        return self.label.format(self.arg)

    def __init__(self, arg):
        self.arg = arg
