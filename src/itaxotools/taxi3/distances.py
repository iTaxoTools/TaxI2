from __future__ import annotations

import re
from pathlib import Path
from typing import NamedTuple
from itertools import chain
from .sequences import Sequence
from .types import Container, Type
from itaxotools.taxi3.library import calculate_distances as calc
from math import isnan, isinf
import alfpy.bbc as bbc
import alfpy.ncd as ncd
from alfpy.utils.seqrecords import SeqRecords


class Distance(NamedTuple):
    metric: DistanceMetric
    x: Sequence
    y: Sequence
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

    def iter_write(self, distances: iter[Distance], *args, **kwargs) -> iter[Distance]:
        raise NotImplementedError()

    def write(self, distances: iter[Distance], *args, **kwargs) -> None:
        for _ in self.iter_write(distances, *args, **kwargs):
            pass


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
            metrics = [DistanceMetric.fromLabel(label) for label in labels]
            for line in f:
                lineData = line[:-1].split('\t')
                idx, idy, labelDistances = lineData[0], lineData[1], lineData[2:]
                distances = (self.distanceFromText(d) for d in labelDistances)
                for distance, metric in zip(distances, metrics):
                    yield Distance(metric, Sequence(idx, None), Sequence(idy, None), distance)

    def iter_write(self, distances: iter[Distance], *args, **kwargs) -> iter[Distance]:
        with open(self.path, 'w') as f:
            metrics = []
            buffer = []
            for d in distances:
                buffer.append(d)
                if len(buffer) > 2:
                    if (buffer[-2].x.id, buffer[-2].y.id) != (d.x.id, d.y.id):
                        break

            for d in buffer[:-1]:
                if str(d.metric) not in metrics:
                    metrics.append(str(d.metric))

            metricString = '\t'.join(metrics)
            metricCount = len(metrics)
            f.write(f'idx\tidy\t{metricString}\n')
            scores = []
            for distance in chain(buffer, distances):
                scores.append(str(distance.d) if distance.d is not None else self.MISSING)
                if len(scores) == metricCount:
                    score = '\t'.join(scores)
                    f.write(f'{distance.x.id}\t{distance.y.id}\t{score}\n')
                    scores = []
                yield distance


class Matrix(DistanceFile):
    MISSING = 'NA'

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
                    yield Distance(metric, Sequence(idx, None), Sequence(idy, None), self.distanceFromText(label_distance))

    def iter_write(self, distances: iter[Distance], *args, **kwargs) -> iter[Distance]:
        with open(self.path, 'w') as f:
            buffer = []

            id = {'idx': [], 'idy': []}
            for distance in distances:
                buffer.append(distance)
                if len(buffer) > 2:
                    if (buffer[-2].x.id) != (distance.x.id):
                        break
                if distance.x.id not in id['idx']:
                    id['idx'].append(distance.x.id)
                if distance.y.id not in id['idy']:
                    id['idy'].append(distance.y.id)

            idy_header = '\t'.join(id['idy'])
            f.write(f'\t{idy_header}\n')
            scores = []
            for distance in chain(buffer, distances):
                d = str(distance.d) if distance.d is not None else self.MISSING
                scores.append(d)
                if len(scores) == len(id['idy']):
                    count = len(id['idy'])
                    score = '\t'.join(scores)
                    f.write(f'{distance.x.id}\t{score}\n')
                    scores = []
                yield distance


class LinearWithExtras(DistanceFile):
    MISSING = 'NA'

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

    def read(
        self,
        idxHeader: str = None,
        idyHeader: str = None,
        tagX: str = ' (query)',
        tagY: str = ' (reference)',
        idxColumn: int = 0,
        idyColumn: int = 1,
    ) -> iter[Distance]:
        with open(self.path, 'r') as f:
            data = f.readline()
            indexLabelStart = 0
            dataList = data.strip().split('\t')
            # print(dataList)
            if idxHeader and idyHeader:
                idxHeader = idxHeader + tagX
                idyHeader = idyHeader + tagY
                idxColumn, idyColumn = dataList.index(idxHeader), dataList.index(idyHeader)

            for label in dataList:
                if not DistanceMetric.fromLabel(label):
                    indexLabelStart +=1
                else:
                    break

            labels = dataList[indexLabelStart:]
            metrics = [DistanceMetric.fromLabel(label) for label in labels]
            extrasHeaderX, extrasHeaderY = dataList[idxColumn+1:idyColumn], dataList[idyColumn+1:indexLabelStart]

            extrasHeaderX = [tag.removesuffix(tagX) for tag in extrasHeaderX]
            extrasHeaderY = [tag.removesuffix(tagY) for tag in extrasHeaderY]

            for line in f:
                lineData = line[:-1].split('\t')
                idx, idy, labelDistances = lineData[idxColumn], lineData[idyColumn], lineData[indexLabelStart:]
                extraDataX = lineData[idxColumn+1:idyColumn]
                extraDataY = lineData[idyColumn+1:indexLabelStart]
                distances = (self.distanceFromText(d) for d in labelDistances)
                extrasX = {k: v for k, v in zip(extrasHeaderX, extraDataX)}
                extrasY = {k: v for k, v in zip(extrasHeaderY, extraDataY)}
                for distance, metric in zip(distances, metrics):
                    yield Distance(metric, Sequence(idx, None, extrasX), Sequence(idy, None, extrasY), distance)

    def iter_write(
        self,
        distances: iter[Distance],
        idxHeader: str = 'seqid',
        idyHeader: str = 'seqid',
        tagX: str = ' (query)',
        tagY: str = ' (reference)',
    ) -> iter[Distance]:
        with open(self.path, 'w') as f:
            metrics = []
            buffer = []

            for d in distances:
                buffer.append(d)
                if len(buffer) > 2:
                    if (buffer[-2].x.id, buffer[-2].y.id) != (d.x.id, d.y.id):
                        break

            for d in buffer[:-1]:
                if str(d.metric) not in metrics:
                    metrics.append(str(d.metric))
            metricCount = len(metrics)

            extrasHeaderX = buffer[0].x.extras.keys()
            extrasHeaderY = buffer[0].y.extras.keys()
            headersX = [idxHeader] + list(extrasHeaderX)
            headersY = [idyHeader] + list(extrasHeaderY)
            headersX = [x + tagX for x in headersX]
            headersY = [y + tagY for y in headersY]
            headers = headersX + headersY + metrics
            headerString = '\t'.join(headers)

            f.write(f'{headerString}\n')

            scores = []
            for distance in chain(buffer, distances):
                scores.append(str(distance.d) if distance.d is not None else self.MISSING)
                if len(scores) == metricCount:
                    score = '\t'.join(scores)
                    extrasX = '\t'.join([distance.x.extras[k] for k in extrasHeaderX])
                    extrasY = '\t'.join([distance.y.extras[k] for k in extrasHeaderY])
                    f.write(f'{distance.x.id}\t{extrasX}\t{distance.y.id}\t{extrasY}\t{score}\n')
                    scores = []
                yield distance


class DistanceMetric(Type):
    """Metrics for calculating distances"""
    label: str

    def __str__(self):
        return self.label

    @staticmethod
    def _is_number(x):
        return not (x is None or isnan(x) or isinf(x))

    def _calculate(self, x: str, y: str) -> float:
        raise NotImplementedError()

    def calculate(self, x: Sequence, y: Sequence) -> Distance:
        return Distance(self, x, y, self._calculate(x.seq, y.seq))

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

    def _calculate(self, x: str, y: str) -> float:
        distance = calc.seq_distances_p(x, y)
        return distance if self._is_number(distance) else None


class UncorrectedWithGaps(DistanceMetric):
    label = 'p-distance with gaps'

    def _calculate(self, x: str, y: str) -> float:
        distance = calc.seq_distances_p_gaps(x, y)
        return distance if self._is_number(distance) else None


class JukesCantor(DistanceMetric):
    label = 'jc'

    def _calculate(self, x: str, y: str) -> float:
        distance = calc.seq_distances_jukes_cantor(x, y)
        return distance if self._is_number(distance) else None


class Kimura2P(DistanceMetric):
    label = 'k2p'

    def _calculate(self, x: str, y: str) -> float:
        distance = calc.seq_distances_kimura2p(x, y)
        return distance if self._is_number(distance) else None


class NCD(DistanceMetric):
    label = 'ncd'

    def _calculate(self, x: str, y: str) -> float:
        records = SeqRecords((0, 1), (x, y))
        distance = ncd.Distance(records)
        d = distance.pairwise_distance(0, 1)
        return d if self._is_number(d) else None


class BBC(DistanceMetric):
    label = 'bbc({})'

    def __init__(self, k=10):
        self.k = k

    def __str__(self):
        return self.label.format(self.k)

    def __eq__(self, other):
        return super().__eq__(other) and self.k == other.k

    def _calculate(self, x: str, y: str) -> float:
        records = SeqRecords((0, 1), (x, y))
        try:
            vector = bbc.create_vectors(records, k=self.k)
            distance = bbc.Distance(vector)
            d = distance.pairwise_distance(0, 1)
        except Exception:
            d = None
        return d if self._is_number(d) else None
