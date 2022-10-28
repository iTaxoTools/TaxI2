from __future__ import annotations

import re
from pathlib import Path
from typing import NamedTuple

from .sequences import Sequence
from .types import Container, Type
from itaxotools.taxi3.library import calculate_distances as calc


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
    MISSING = 'NA'
    dataList = []

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
                    yield Distance(metric, idx, idy, self.distanceFromText(labelDistances[label]))

    def write(self, distances: iter[Distance], *args, **kwargs) -> None:
        metrics = []
        # get metrics
        self.dataList = []
        for d in distances:
            if str(d.metric) not in metrics:
                metrics.append(str(d.metric))
            self.dataList.append(d)
        with open(self.path, 'w') as f:
            metricString = '\t'.join(metrics)
            metricCount = len(metrics)
            f.write(f'idx\tidy\t{metricString}\n')
            count = metricCount
            print('ewww')
            if metricCount > 1:
                scores = []
                for distance in self.dataList:
                    print(d)
                    scores.append(str(distance.d) if distance.d else self.MISSING)
                    count -= 1
                    if len(scores) == metricCount:
                        count = metricCount
                        score = '\t'.join(scores)
                        print(f'{distance.idx}\t{distance.idy}\t{score}\n')
                        f.write(f'{distance.idx}\t{distance.idy}\t{score}\n')
                        scores = []
            else:
                for distance in self.dataList:
                    print('em')
                    score = str(distance.d) if distance.d is not None else self.MISSING
                    f.write(f'{distance.idx}\t{distance.idy}\t{score}\n')
            f.close()


class Matrix(DistanceFile):
    MISSING = 'NA'
    dataList = []

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
            for line in f:
                data = line.strip().split('\t')
                idx = data[0]
                for index in range(len(data[1:])):
                    idy = id_to_compare[index]
                    label_distance = data[1:][index]
                    yield Distance(metric, idx, idy, self.distanceFromText(label_distance))

    def write(self, distances: iter[Distance], *args, **kwargs) -> None:
        self.dataList = []
        id = {'idx': [], 'idy': []}
        for distance in distances:
            if distance.idx not in id['idx']:
                id['idx'].append(distance.idx)
            if distance.idy not in id['idy']:
                id['idy'].append(distance.idy)
            self.dataList.append(distance)

        with open(self.path, 'w') as f:
            idy_header = '\t'.join(id['idy'])
            f.write(f'\t{idy_header}\n')
            count = len(id['idy'])
            scores = []
            for distance in self.dataList:
                d = str(distance.d) if distance.d is not None else self.MISSING
                scores.append(d)
                count -= 1
                if len(scores) == len(id['idy']):
                    count = len(id['idy'])
                    score = '\t'.join(scores)
                    f.write(f'{distance.idx}\t{score}\n')
                    scores = []


class DistanceMetric(Type):
    """Metrics for calculating distances"""
    label: str

    def __str__(self):
        return self.label

    def _calculate(self, x: str, y: str) -> float:
        distance = None
        if self.label == 'p-distance':
            distance = calc.seq_distances_p(x, y)
        elif self.label == 'p-distance with gaps':
            print('entered')
            distance = calc.seq_distances_p_gaps(x, y)
        elif self.label == 'jc':
            distance = calc.seq_distances_jukes_cantor(x, y)
        elif self.label == 'k2p':
            distance = calc.seq_distances_kimura2p(x, y)

        return distance if str(distance) != 'nan' else None

    def calculate(self, x: Sequence, y: Sequence) -> Distance:
        return Distance(self, x.id, y.id, self._calculate(x.seq, y.seq))

    @classmethod
    def fromLabel(cls, label: str):
        label_arg = None
        res = re.search(r'(\w+)\((\d+)\)', label)
        if res:
            label = res.group(1) + '({})'
            label_arg = res.group(2)
        for child in cls:
            if label == child.label:
                if label_arg:
                    return child(int(label_arg))
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
    label = 'ncd'


class BBC(DistanceMetric):
    label = 'bbc({})'

    def __init__(self, arg):
        self.arg = arg

    def __str__(self):
        return self.label.format(self.arg)

    def __eq__(self, other):
        return super().__eq__(other) and self.arg == other.arg
