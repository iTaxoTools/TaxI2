from __future__ import annotations

from pathlib import Path
from sys import stderr
from typing import Callable, NamedTuple

import pytest
from utility import assert_eq_files

from itaxotools.taxi3.distances import (
    Distance, DistanceFile, DistanceMetric, Distances)
from itaxotools.taxi3.sequences import Sequence

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class ReadTest(NamedTuple):
    fixture: Callable[[], Distances]
    input: str
    file: DistanceFile
    kwargs: dict = {}

    def validate(self, generated: Distances):
        fixture_list = list(self.fixture())
        generated_list = list(generated)
        for distance in fixture_list:
            assert distance in generated_list


class WriteTest(NamedTuple):
    fixture: Callable[[], Distances]
    output: str
    file: DistanceFile
    kwargs: dict = {}

    def generate(self) -> Distances:
        return self.fixture()


class LabelTest(NamedTuple):
    metric: DistanceMetric
    label: str

    def check(self):
        assert self.metric == DistanceMetric.fromLabel(self.label)
        assert self.label == str(self.metric)


class MetricTest(NamedTuple):
    metric: DistanceMetric
    seq_x: str
    seq_y: str
    d: float
    precision: float = 0.0

    def check(self):
        x = Sequence('idx', self.seq_x)
        y = Sequence('idy', self.seq_y)
        r = self.metric.calculate(x, y)
        assert r.metric == self.metric
        assert r.x.id == 'idx'
        assert r.y.id == 'idy'
        if isinstance(r.d, float):
            assert abs(r.d - self.d) <= self.precision
        else:
            assert r.d == self.d
            assert r.d is None


class MetricFileTest(NamedTuple):
    file: str
    precision: float

    def get_metric_tests(self) -> iter[MetricTest]:
        path = TEST_DATA_DIR / self.file
        with open(path, 'r') as f:
            data = f.readline()
            labels = data.strip().split('\t')[2:]
            metrics = [DistanceMetric.fromLabel(label) for label in labels]
            for line in f:
                lineData = line[:-1].split('\t')
                seqX, seqY, labelDistances = lineData[0], lineData[1], lineData[2:]
                distances = (DistanceFile.Linear.distanceFromText(d) for d in labelDistances)
                for distance, metric in zip(distances, metrics):
                    yield MetricTest(metric, seqX, seqY, distance, self.precision)


def distances_simple() -> Distances:
    metric = DistanceMetric.Uncorrected()
    return Distances([
        Distance(metric, Sequence('id1', None), Sequence('id2', None), 0.1),
        Distance(metric, Sequence('id1', None), Sequence('id3', None), 0.2),
        Distance(metric, Sequence('id1', None), Sequence('id4', None), 0.3),
    ])


def distances_multiple() -> Distances:
    return Distances([
        Distance(DistanceMetric.Uncorrected(), Sequence('id1', None), Sequence('id2', None), 0.11),
        Distance(DistanceMetric.UncorrectedWithGaps(), Sequence('id1', None), Sequence('id2', None), 0.12),
        Distance(DistanceMetric.JukesCantor(), Sequence('id1', None), Sequence('id2', None), 0.13),
        Distance(DistanceMetric.Kimura2P(), Sequence('id1', None), Sequence('id2', None), 0.14),
        Distance(DistanceMetric.NCD(), Sequence('id1', None), Sequence('id2', None), 0.15),
        Distance(DistanceMetric.BBC(0), Sequence('id1', None), Sequence('id2', None), 0.16),

        Distance(DistanceMetric.Uncorrected(), Sequence('id1', None), Sequence('id3', None), 0.21),
        Distance(DistanceMetric.UncorrectedWithGaps(), Sequence('id1', None), Sequence('id3', None), 0.22),
        Distance(DistanceMetric.JukesCantor(), Sequence('id1', None), Sequence('id3', None), 0.23),
        Distance(DistanceMetric.Kimura2P(), Sequence('id1', None), Sequence('id3', None), 0.24),
        Distance(DistanceMetric.NCD(), Sequence('id1', None), Sequence('id3', None), 0.25),
        Distance(DistanceMetric.BBC(0), Sequence('id1', None), Sequence('id3', None), 0.26),

        Distance(DistanceMetric.Uncorrected(), Sequence('id1', None), Sequence('id4', None), 0.31),
        Distance(DistanceMetric.UncorrectedWithGaps(), Sequence('id1', None), Sequence('id4', None), 0.32),
        Distance(DistanceMetric.JukesCantor(), Sequence('id1', None), Sequence('id4', None), 0.33),
        Distance(DistanceMetric.Kimura2P(), Sequence('id1', None), Sequence('id4', None), 0.34),
        Distance(DistanceMetric.NCD(), Sequence('id1', None), Sequence('id4', None), 0.35),
        Distance(DistanceMetric.BBC(0), Sequence('id1', None), Sequence('id4', None), 0.36),
    ])


def distances_square() -> Distances:
    metric = DistanceMetric.Uncorrected()
    return Distances([
        Distance(metric, Sequence('id1', None), Sequence('id1', None), 0.0),
        Distance(metric, Sequence('id1', None), Sequence('id2', None), 0.1),
        Distance(metric, Sequence('id1', None), Sequence('id3', None), 0.2),

        Distance(metric, Sequence('id2', None), Sequence('id1', None), 0.1),
        Distance(metric, Sequence('id2', None), Sequence('id2', None), 0.0),
        Distance(metric, Sequence('id2', None), Sequence('id3', None), 0.3),

        Distance(metric, Sequence('id3', None), Sequence('id1', None), 0.2),
        Distance(metric, Sequence('id3', None), Sequence('id2', None), 0.3),
        Distance(metric, Sequence('id3', None), Sequence('id3', None), 0.0),
    ])


def distances_square_unknown() -> Distances:
    metric = DistanceMetric.Unknown()
    return Distances([
        Distance(metric, dis.x, dis.y, dis.d) for dis in distances_square()
    ])


def distances_rectangle() -> Distances:
    metric = DistanceMetric.Uncorrected()
    return Distances([
        Distance(metric, Sequence('id1', None), Sequence('id4', None), 0.14),
        Distance(metric, Sequence('id1', None), Sequence('id5', None), 0.15),
        Distance(metric, Sequence('id1', None), Sequence('id6', None), 0.16),
        Distance(metric, Sequence('id1', None), Sequence('id7', None), 0.17),
        Distance(metric, Sequence('id1', None), Sequence('id8', None), 0.18),
        Distance(metric, Sequence('id1', None), Sequence('id9', None), 0.19),

        Distance(metric, Sequence('id2', None), Sequence('id4', None), 0.24),
        Distance(metric, Sequence('id2', None), Sequence('id5', None), 0.25),
        Distance(metric, Sequence('id2', None), Sequence('id6', None), 0.26),
        Distance(metric, Sequence('id2', None), Sequence('id7', None), 0.27),
        Distance(metric, Sequence('id2', None), Sequence('id8', None), 0.28),
        Distance(metric, Sequence('id2', None), Sequence('id9', None), 0.29),

        Distance(metric, Sequence('id3', None), Sequence('id4', None), 0.34),
        Distance(metric, Sequence('id3', None), Sequence('id5', None), 0.35),
        Distance(metric, Sequence('id3', None), Sequence('id6', None), 0.36),
        Distance(metric, Sequence('id3', None), Sequence('id7', None), 0.37),
        Distance(metric, Sequence('id3', None), Sequence('id8', None), 0.38),
        Distance(metric, Sequence('id3', None), Sequence('id9', None), 0.39),
    ])


def distances_missing() -> Distances:
    metric = DistanceMetric.Uncorrected()
    return Distances([
        Distance(metric, Sequence('id1', None), Sequence('id1', None), 0.0),
        Distance(metric, Sequence('id1', None), Sequence('id2', None), None),

        Distance(metric, Sequence('id2', None), Sequence('id1', None), None),
        Distance(metric, Sequence('id2', None), Sequence('id2', None), 0.0),
    ])


def distances_extras() -> Distances:
    return Distances([
        Distance(
            DistanceMetric.Uncorrected(),
            Sequence('query1', None, dict(voucher='K')),
            Sequence('reference1', None, dict(voucher='X', organism='A')),
            0.11),
        Distance(
            DistanceMetric.UncorrectedWithGaps(),
            Sequence('query1', None, dict(voucher='K')),
            Sequence('reference1', None, dict(voucher='X', organism='A')),
            0.12),
        Distance(
            DistanceMetric.JukesCantor(),
            Sequence('query1', None, dict(voucher='K')),
            Sequence('reference1', None, dict(voucher='X', organism='A')),
            0.13),
        Distance(
            DistanceMetric.Kimura2P(),
            Sequence('query1', None, dict(voucher='K')),
            Sequence('reference1', None, dict(voucher='X', organism='A')),
            0.14),

        Distance(
            DistanceMetric.Uncorrected(),
            Sequence('query1', None, dict(voucher='K')),
            Sequence('reference2', None, dict(voucher='Y', organism='B')),
            0.21),
        Distance(
            DistanceMetric.UncorrectedWithGaps(),
            Sequence('query1', None, dict(voucher='K')),
            Sequence('reference2', None, dict(voucher='Y', organism='B')),
            0.22),
        Distance(
            DistanceMetric.JukesCantor(),
            Sequence('query1', None, dict(voucher='K')),
            Sequence('reference2', None, dict(voucher='Y', organism='B')),
            0.23),
        Distance(
            DistanceMetric.Kimura2P(),
            Sequence('query1', None, dict(voucher='K')),
            Sequence('reference2', None, dict(voucher='Y', organism='B')),
            0.24),

        Distance(
            DistanceMetric.Uncorrected(),
            Sequence('query2', None, dict(voucher='L')),
            Sequence('reference3', None, dict(voucher='Z', organism='C')),
            0.31),
        Distance(
            DistanceMetric.UncorrectedWithGaps(),
            Sequence('query2', None, dict(voucher='L')),
            Sequence('reference3', None, dict(voucher='Z', organism='C')),
            0.32),
        Distance(
            DistanceMetric.JukesCantor(),
            Sequence('query2', None, dict(voucher='L')),
            Sequence('reference3', None, dict(voucher='Z', organism='C')),
            0.33),
        Distance(
            DistanceMetric.Kimura2P(),
            Sequence('query2', None, dict(voucher='L')),
            Sequence('reference3', None, dict(voucher='Z', organism='C')),
            None),
    ])


read_tests = [
    ReadTest(distances_simple, 'simple.linear', DistanceFile.Linear),
    ReadTest(distances_multiple, 'multiple.linear', DistanceFile.Linear),
    ReadTest(distances_missing, 'missing.linear', DistanceFile.Linear),

    ReadTest(distances_square_unknown, 'square.matrix', DistanceFile.Matrix),
    ReadTest(distances_square, 'square.matrix', DistanceFile.Matrix,
        dict(metric=DistanceMetric.Uncorrected())),
    ReadTest(distances_rectangle, 'rectangle.matrix', DistanceFile.Matrix,
        dict(metric=DistanceMetric.Uncorrected())),
    ReadTest(distances_missing, 'missing.matrix', DistanceFile.Matrix,
        dict(metric=DistanceMetric.Uncorrected())),

    ReadTest(distances_extras, 'extras.tsv', DistanceFile.LinearWithExtras,
        dict(idxHeader='seqid', idyHeader='id', tagX='_x', tagY='_y')),
    ReadTest(distances_extras, 'extras.tsv', DistanceFile.LinearWithExtras,
        dict(idxColumn=0, idyColumn=2, tagX='_x', tagY='_y')),
]


write_tests = [
    WriteTest(distances_simple, 'simple.linear', DistanceFile.Linear),
    WriteTest(distances_multiple, 'multiple.linear', DistanceFile.Linear),
    WriteTest(distances_missing, 'missing.linear', DistanceFile.Linear),
    WriteTest(distances_square, 'square.matrix', DistanceFile.Matrix),
    WriteTest(distances_rectangle, 'rectangle.matrix', DistanceFile.Matrix),
    WriteTest(distances_missing, 'missing.matrix', DistanceFile.Matrix),
    WriteTest(distances_extras, 'extras.tsv', DistanceFile.LinearWithExtras,
        dict(idxHeader='seqid', idyHeader='id', tagX='_x', tagY='_y')),
    WriteTest(distances_missing, 'missing.formatted.linear', DistanceFile.Linear,
        dict(format='{:.2e}', missing='nan')),
    WriteTest(distances_missing, 'missing.formatted.linear', DistanceFile.LinearWithExtras,
        dict(format='{:.2e}', missing='nan')),
    WriteTest(distances_missing, 'missing.formatted.matrix', DistanceFile.Matrix,
        dict(format='{:.2e}', missing='nan')),
]


label_tests = [
    LabelTest(DistanceMetric.Uncorrected(), 'p-distance'),
    LabelTest(DistanceMetric.UncorrectedWithGaps(), 'p-distance with gaps'),
    LabelTest(DistanceMetric.JukesCantor(), 'jc'),
    LabelTest(DistanceMetric.Kimura2P(), 'k2p'),
    LabelTest(DistanceMetric.NCD(), 'ncd'),
    LabelTest(DistanceMetric.NCD(), 'ncd'),
    LabelTest(DistanceMetric.BBC(0), 'bbc(0)'),
    LabelTest(DistanceMetric.BBC(1), 'bbc(1)'),
]


metric_tests = [
    MetricTest(DistanceMetric.Uncorrected(), 'gg-ccnccta', 'ggaccaccaa', 1.0 / 8.0),
    MetricTest(DistanceMetric.UncorrectedWithGaps(), 'gg-ccnccta', 'ggaccaccaa', 2.0 / 9.0),
    MetricTest(DistanceMetric.Uncorrected(), '---', 'nnn', None),
]


metric_file_tests = [
    MetricFileTest('metrics.tsv', 0.00051),
]


@pytest.mark.parametrize("test", read_tests)
def test_read_distances(test: ReadTest) -> None:
    input_path = TEST_DATA_DIR / test.input
    distances = test.file(input_path).read(**test.kwargs)
    test.validate(distances)


@pytest.mark.parametrize("test", write_tests)
def test_write_distances(test: WriteTest, tmp_path: Path) -> None:
    fixed_path = TEST_DATA_DIR / test.output
    output_path = tmp_path / test.output
    distances = iter(test.generate())
    test.file(output_path).write(distances, **test.kwargs)
    assert_eq_files(output_path, fixed_path)


@pytest.mark.parametrize("test", label_tests)
def test_labels(test: LabelTest) -> None:
    test.check()


@pytest.mark.parametrize("test", metric_tests)
def test_metrics(test: MetricTest) -> None:
    test.check()


@pytest.mark.parametrize("test", metric_file_tests)
def test_metrics_from_files(test: MetricFileTest) -> None:
    stack = []
    for metric_test in test.get_metric_tests():
        try:
            metric_test.check()
        except AssertionError as a:
            stack.append(a)
    for a in stack:
        print(a.args[0], '\n', file=stderr)
    assert len(stack) == 0
