from __future__ import annotations

import re
from itertools import chain
from math import isinf, isnan
from pathlib import Path
from typing import NamedTuple, Generator

import alfpy.bbc as bbc
import alfpy.ncd as ncd
from alfpy.utils.seqrecords import SeqRecords

from itaxotools.taxi3.library import calculate_distances as calc

from .sequences import Sequence
from .types import Container, Type
from .handlers import FileHandler, ReadHandle, WriteHandle


class Distance(NamedTuple):
    metric: DistanceMetric
    x: Sequence
    y: Sequence
    d: float | None


class Distances(Container[Distance]):
    @classmethod
    def fromPath(cls, path: Path, handler: DistanceHandler, *args, **kwargs) -> Distances:
        return cls(handler, path, *args, **kwargs)


class DistanceHandler(FileHandler[Distance]):
    def _open(
        self,
        path: Path,
        mode: 'r' | 'w' = 'r',
        missing: str = 'NA',
        formatter: str = '{:f}',
        *args, **kwargs
    ):
        self.missing = missing
        self.formatter = formatter
        super()._open(path, mode, *args, **kwargs)

    def distanceFromText(self, text: str) -> float | None:
        if text == self.missing:
            return None
        return float(text)

    def distanceToText(self, d: float | None) -> str:
        if d is None:
            return self.missing
        return self.formatter.format(d)


class Linear(DistanceHandler):
    def _open_writable(self, *args, **kwargs):
        self.buffer: list[Distance] = []
        self.wrote_headers = False
        super()._open_writable()

    def _iter_read(self) -> ReadHandle[Distance]:
        with FileHandler.Tabfile(self.path, 'r', has_headers=True) as file:
            metrics = [DistanceMetric.fromLabel(label) for label in file.headers[2:]]
            yield self
            for row in file:
                idx, idy, distances = row[0], row[1], row[2:]
                distances = (self.distanceFromText(d) for d in distances)
                for distance, metric in zip(distances, metrics):
                    yield Distance(metric, Sequence(idx, None), Sequence(idy, None), distance)

    def _iter_write(self) -> WriteHandle[Distance]:
        with FileHandler.Tabfile(self.path, 'w') as file:
            try:
                line = yield from self._assemble_line()
                self._write_headers(file, line)
                self._write_scores(file, line)
                while True:
                    line = yield from self._assemble_line()
                    self._write_scores(file, line)
            except GeneratorExit:
                line = self.buffer
                self._write_headers(file, line)
                self._write_scores(file, line)
                return

    def _assemble_line(self) -> Generator[None, Distance, list[Distance]]:
        buffer = self.buffer
        try:
            while True:
                distance = yield
                buffer.append(distance)
                if any((
                    buffer[0].x.id != buffer[-1].x.id,
                    buffer[0].y.id != buffer[-1].y.id,
                )):
                    self.buffer = buffer[-1:]
                    return buffer[:-1]
        except GeneratorExit:
            return

    def _write_headers(self, file: FileHandler.Tabfile, line: list[Distance]):
        if self.wrote_headers:
            return
        metrics = [str(distance.metric) for distance in line]
        out = ('idx', 'idy', *metrics)
        file.write(out)
        self.wrote_headers = True

    def _write_scores(self, file: FileHandler.Tabfile, line: list[Distance]):
        scores = [self.distanceToText(distance.d) for distance in line]
        out = (line[0].x.id, line[0].y.id, *scores)
        file.write(out)


class Matrix(DistanceHandler):
    def _open_readable(self, metric: DistanceMetric = None, *args, **kwargs):
        self.metric = metric or DistanceMetric.Unknown()
        super()._open_readable()

    def _open_writable(self, *args, **kwargs):
        self.buffer: list[Distance] = []
        self.wrote_headers = False
        super()._open_writable()

    def _iter_read(self) -> ReadHandle[Distance]:
        with FileHandler.Tabfile(self.path, 'r', has_headers=True) as file:
            idys = file.headers[1:]
            metric = self.metric
            yield self
            for row in file:
                idx, scores = row[0], row[1:]
                seqx = Sequence(idx, None)
                for score, idy in zip(scores, idys):
                    d = self.distanceFromText(score)
                    yield Distance(metric, seqx, Sequence(idy, None), d)

    def _iter_write(self) -> WriteHandle[Distance]:
        with FileHandler.Tabfile(self.path, 'w') as file:
            try:
                line = yield from self._assemble_line()
                self._write_headers(file, line)
                self._write_scores(file, line)
                while True:
                    line = yield from self._assemble_line()
                    self._write_scores(file, line)
            except GeneratorExit:
                line = self.buffer
                self._write_headers(file, line)
                self._write_scores(file, line)
                return

    def _assemble_line(self) -> Generator[None, Distance, list[Distance]]:
        buffer = self.buffer
        try:
            while True:
                distance = yield
                buffer.append(distance)
                if buffer[0].x.id != buffer[-1].x.id:
                    self.buffer = buffer[-1:]
                    return buffer[:-1]
        except GeneratorExit:
            return

    def _write_headers(self, file: FileHandler.Tabfile, line: list[Distance]):
        if self.wrote_headers:
            return
        idys = [distance.y.id for distance in line]
        out = ('', *idys)
        file.write(out)
        self.wrote_headers = True

    def _write_scores(self, file: FileHandler.Tabfile, line: list[Distance]):
        scores = [self.distanceToText(distance.d) for distance in line]
        out = (line[0].x.id, *scores)
        file.write(out)


class LinearWithExtras(DistanceHandler):
    # pending rewrite
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
                    indexLabelStart += 1
                else:
                    break

            labels = dataList[indexLabelStart:]
            metrics = [DistanceMetric.fromLabel(label) for label in labels]
            extrasHeaderX, extrasHeaderY = dataList[idxColumn + 1:idyColumn], dataList[idyColumn + 1:indexLabelStart]

            extrasHeaderX = [tag.removesuffix(tagX) for tag in extrasHeaderX]
            extrasHeaderY = [tag.removesuffix(tagY) for tag in extrasHeaderY]

            for line in f:
                lineData = line[:-1].split('\t')
                idx, idy, labelDistances = lineData[idxColumn], lineData[idyColumn], lineData[indexLabelStart:]
                extraDataX = lineData[idxColumn + 1:idyColumn]
                extraDataY = lineData[idyColumn + 1:indexLabelStart]
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
        *args,
        **kwargs) -> iter[Distance]:
        with open(self.path, 'w') as f:
            metrics = []
            buffer = []
            self.MISSING = kwargs.get('missing', 'NA')
            scroreFormat = kwargs.get('formatScore', None)

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
                scores.append(str(scroreFormat.format(float(distance.d)) if scroreFormat else distance.d) if distance.d is not None else self.MISSING)
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
