from __future__ import annotations

import re
from typing import Callable, Iterable, NamedTuple
from enum import Enum, auto
import numpy as np
from .sequences import Sequence, Sequences
from .types import Type, Container
import parse

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
    MISSING = 'nan'

    @classmethod
    def distanceFromText(cls, text: str) -> float | None:
        if text == cls.MISSING:
            return None
        return float(text)

    @classmethod
    def distanceToText(cls, d: float | None) -> str:
        if d is None:
            return cls.MISSING
        return str(d)

    def read(self) -> iter[Distance]:
        with open(self.path, 'r') as f:
            data = f.readline()
            labels = data.strip().split('\t')[2:]
            for line in f:
                lineData = line[:-1].split('\t')
                idx, idy, labelDistances = lineData[0], lineData[1], lineData[2:]
                for label in range(len(labels)):
                    metric = DistanceMetric.fromLabel(labels[label])
                    #print(Distance(metric, idx, idy, self.distanceFromText(labelDistances[label])))
                    yield Distance(metric, idx, idy, self.distanceFromText(labelDistances[label]))


class Matrix(DistanceFile):
    MISSING = 'nan'

    @classmethod
    def distanceFromText(cls, text: str) -> float | None:
        if text == cls.MISSING:
            return None
        return float(text)

    def read(self, metric: DistanceMetric = None) -> iter[Distance]:
        metric = metric or DistanceMetric.Unknown()
        with open(self.path, 'r') as f:
            id_data = f.readline()
            id_to_compare = id_data.strip().split('\t')
            print(id_to_compare)
            for line in f:
                data = line.strip().split('\t')
                idx = data[0]
                for index in range(len(data[1:])):
                    idy = id_to_compare[index]
                    label_distance = data[1:][index]
                    yield Distance(metric, idx, idy, self.distanceFromText(label_distance))


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
        label_arg = None
        res = re.search(r'(\w+)\((\d+)\)', label)
        if res:
            label = res.group(1)+'({})'
            label_arg = res.group(2)
        for child in cls:
            if label == child.label:
                if label_arg:
                    return child(label_arg)
                else:
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

    # def __eq__(self, other):
    #     return super().__eq__(other) and self.arg == other.arg


class BBC(DistanceMetric):
    label = 'bbc({})'

    def __init__(self, arg):
        self.arg = arg

    def __str__(self):
        return self.label.format(self.arg)

    # def __eq__(self, other):
    #     return super().__eq__(other) and self.arg == other.arg
