from __future__ import annotations

import re
from math import isinf, isnan
from pathlib import Path
from typing import Generator, Literal, NamedTuple

import alfpy.bbc as bbc
import alfpy.ncd as ncd
from alfpy.utils.seqrecords import SeqRecords

from itaxotools import calculate_distances as calc

from .handlers import FileHandler, ReadHandle, WriteHandle
from .sequences import Sequence
from .types import Container, Type


class Distance(NamedTuple):
    metric: DistanceMetric
    x: Sequence
    y: Sequence
    d: float | None


class Distances(Container[Distance]):
    @classmethod
    def fromPath(
        cls, path: Path, handler: DistanceHandler, *args, **kwargs
    ) -> Distances:
        return cls(handler, path, *args, **kwargs)


class DistanceHandler(FileHandler[Distance]):
    def _open(
        self,
        path: Path,
        mode: Literal["r", "w"] = "r",
        missing: str = "NA",
        formatter: str = "{:f}",
        *args,
        **kwargs,
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
    def _iter_read(self) -> ReadHandle[Distance]:
        with FileHandler.Tabfile(self.path, "r", has_headers=True) as file:
            if file.headers is None:
                yield self
                return
            metrics = [DistanceMetric.fromLabel(label) for label in file.headers[2:]]
            yield self
            for row in file:
                idx, idy, distances = row[0], row[1], row[2:]
                distances = (self.distanceFromText(d) for d in distances)
                for distance, metric in zip(distances, metrics):
                    yield Distance(
                        metric, Sequence(idx, None), Sequence(idy, None), distance
                    )

    def _iter_write(self) -> WriteHandle[Distance]:
        self.buffer: list[Distance] = []
        self.wrote_headers = False

        with FileHandler.Tabfile(self.path, "w") as file:
            try:
                line = yield from self._assemble_line()
                self._write_headers(file, line)
                self._write_scores(file, line)
                while True:
                    line = yield from self._assemble_line()
                    self._write_scores(file, line)
            except GeneratorExit:
                line = self.buffer
                if not line:
                    return
                self._write_headers(file, line)
                self._write_scores(file, line)
                return

    def _assemble_line(self) -> Generator[None, Distance, list[Distance]]:
        buffer = self.buffer
        try:
            while True:
                distance = yield
                buffer.append(distance)
                if any(
                    (
                        buffer[0].x.id != buffer[-1].x.id,
                        buffer[0].y.id != buffer[-1].y.id,
                    )
                ):
                    self.buffer = buffer[-1:]
                    return buffer[:-1]
        except GeneratorExit:
            return

    def _write_headers(self, file: FileHandler.Tabfile, line: list[Distance]):
        if self.wrote_headers:
            return
        metrics = [str(distance.metric) for distance in line]
        out = ("idx", "idy", *metrics)
        file.write(out)
        self.wrote_headers = True

    def _write_scores(self, file: FileHandler.Tabfile, line: list[Distance]):
        scores = [self.distanceToText(distance.d) for distance in line]
        out = (line[0].x.id, line[0].y.id, *scores)
        file.write(out)


class Matrix(DistanceHandler):
    def _iter_read(self, metric: DistanceMetric = None) -> ReadHandle[Distance]:
        metric = metric or DistanceMetric.Unknown()

        with FileHandler.Tabfile(self.path, "r", has_headers=True) as file:
            if file.headers is None:
                yield self
                return
            idys = file.headers[1:]
            yield self
            for row in file:
                idx, scores = row[0], row[1:]
                seqx = Sequence(idx, None)
                for score, idy in zip(scores, idys):
                    d = self.distanceFromText(score)
                    yield Distance(metric, seqx, Sequence(idy, None), d)

    def _iter_write(self) -> WriteHandle[Distance]:
        self.buffer: list[Distance] = []
        self.wrote_headers = False

        with FileHandler.Tabfile(self.path, "w") as file:
            try:
                line = yield from self._assemble_line()
                self._write_headers(file, line)
                self._write_scores(file, line)
                while True:
                    line = yield from self._assemble_line()
                    self._write_scores(file, line)
            except GeneratorExit:
                line = self.buffer
                if not line:
                    return
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
        out = ("", *idys)
        file.write(out)
        self.wrote_headers = True

    def _write_scores(self, file: FileHandler.Tabfile, line: list[Distance]):
        scores = [self.distanceToText(distance.d) for distance in line]
        out = (line[0].x.id, *scores)
        file.write(out)


class WithExtras(DistanceHandler.Linear):
    def _iter_read(
        self,
        idxHeader: str = None,
        idyHeader: str = None,
        tagX: str = " (query)",
        tagY: str = " (reference)",
        idxColumn: int = 0,
        idyColumn: int = 1,
    ) -> ReadHandle[Distance]:
        with FileHandler.Tabfile(self.path, "r", has_headers=True) as file:
            if file.headers is None:
                yield self
                return
            headers = file.headers

            if idxHeader and idyHeader:
                idxHeader = idxHeader + tagX
                idyHeader = idyHeader + tagY
                idxColumn = headers.index(idxHeader)
                idyColumn = headers.index(idyHeader)

            try:
                metricIndexStart = next(
                    i for i, x in enumerate(headers) if DistanceMetric.fromLabel(x)
                )
            except StopIteration:
                raise Exception("No metrics found in the header line!")

            sliceX = slice(idxColumn + 1, idyColumn)
            sliceY = slice(idyColumn + 1, metricIndexStart)

            metricHeaders = headers[metricIndexStart:]
            metrics = [DistanceMetric.fromLabel(header) for header in metricHeaders]
            extrasHeaderX = [header.removesuffix(tagX) for header in headers[sliceX]]
            extrasHeaderY = [header.removesuffix(tagY) for header in headers[sliceY]]

            yield self

            for row in file:
                idx = row[idxColumn]
                idy = row[idyColumn]
                extraDataX = row[sliceX]
                extraDataY = row[sliceY]
                extrasX = {k: v for k, v in zip(extrasHeaderX, extraDataX)}
                extrasY = {k: v for k, v in zip(extrasHeaderY, extraDataY)}
                distances = (self.distanceFromText(d) for d in row[metricIndexStart:])
                for distance, metric in zip(distances, metrics):
                    yield Distance(
                        metric,
                        Sequence(idx, None, extrasX),
                        Sequence(idy, None, extrasY),
                        distance,
                    )

    def _iter_write(
        self,
        idxHeader: str = "seqid",
        idyHeader: str = "seqid",
        tagX: str = " (query)",
        tagY: str = " (reference)",
    ) -> WriteHandle[Distance]:
        self.idxHeader = idxHeader
        self.idyHeader = idyHeader
        self.tagX = tagX
        self.tagY = tagY

        yield from super()._iter_write()

    def _write_headers(self, file: FileHandler.Tabfile, line: list[Distance]):
        if self.wrote_headers:
            return
        idxHeader = self.idxHeader + self.tagX
        idyHeader = self.idyHeader + self.tagY
        extrasX = [key + self.tagX for key in line[0].x.extras.keys()]
        extrasY = [key + self.tagY for key in line[0].y.extras.keys()]
        metrics = [str(distance.metric) for distance in line]
        out = (idxHeader, *extrasX, idyHeader, *extrasY, *metrics)
        file.write(out)
        self.wrote_headers = True

    def _write_scores(self, file: FileHandler.Tabfile, line: list[Distance]):
        idx = line[0].x.id
        idy = line[0].y.id
        extrasX = line[0].x.extras.values()
        extrasY = line[0].y.extras.values()
        extrasX = [x if x is not None else self.missing for x in extrasX]
        extrasY = [y if y is not None else self.missing for y in extrasY]
        scores = [self.distanceToText(distance.d) for distance in line]
        out = (idx, *extrasX, idy, *extrasY, *scores)
        file.write(out)


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
        res = re.search(r"(\w+)\((\d+)\)", label)
        if res:
            label = res.group(1) + "({})"
            label_arg = res.group(2)
        for child in cls:
            if label == child.label:
                if label_arg:
                    return child(int(label_arg))
                else:
                    return child()


class Unknown(DistanceMetric):
    label = "?"


class Uncorrected(DistanceMetric):
    label = "p"

    def _calculate(self, x: str, y: str) -> float:
        distance = calc.seq_distances_p(x, y)
        return distance if self._is_number(distance) else None


class UncorrectedWithGaps(DistanceMetric):
    label = "p-gaps"

    def _calculate(self, x: str, y: str) -> float:
        distance = calc.seq_distances_p_gaps(x, y)
        return distance if self._is_number(distance) else None


class JukesCantor(DistanceMetric):
    label = "jc"

    def _calculate(self, x: str, y: str) -> float:
        distance = calc.seq_distances_jukes_cantor(x, y)
        return distance if self._is_number(distance) else None


class Kimura2P(DistanceMetric):
    label = "k2p"

    def _calculate(self, x: str, y: str) -> float:
        distance = calc.seq_distances_kimura2p(x, y)
        return distance if self._is_number(distance) else None


class NCD(DistanceMetric):
    label = "ncd"

    def _calculate(self, x: str, y: str) -> float:
        records = SeqRecords((0, 1), (x, y))
        distance = ncd.Distance(records)
        d = distance.pairwise_distance(0, 1)
        return d if self._is_number(d) else None


class BBC(DistanceMetric):
    label = "bbc({})"

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
