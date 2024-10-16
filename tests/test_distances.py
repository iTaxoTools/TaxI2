from __future__ import annotations

from pathlib import Path
from sys import stderr
from typing import Callable, Iterator, NamedTuple

import pytest
from utility import assert_eq_files

from itaxotools.taxi2.distances import (
    Distance,
    DistanceHandler,
    DistanceMetric,
    Distances,
)
from itaxotools.taxi2.sequences import Sequence, Sequences

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class ReadTest(NamedTuple):
    fixture: Callable[[], Distances]
    input: str
    handler: DistanceHandler
    kwargs: dict = {}

    @property
    def input_path(self) -> Path:
        return TEST_DATA_DIR / self.input

    @property
    def fixed(self) -> Distances:
        return self.fixture()

    def validate(self) -> None:
        distances = Distances.fromPath(self.input_path, self.handler, **self.kwargs)
        generated_list = list(distances)
        fixed_list = list(self.fixed)
        assert len(fixed_list) == len(generated_list)
        for distance in fixed_list:
            assert distance in generated_list


class WriteTest(NamedTuple):
    fixture: Callable[[], Distances]
    output: str
    handler: DistanceHandler
    kwargs: dict = {}

    @property
    def fixed_path(self) -> Path:
        return TEST_DATA_DIR / self.output

    @property
    def fixed(self) -> Sequences:
        return self.fixture()

    def get_output_path(self, tmp_path) -> Path:
        return tmp_path / self.output

    def validate(self, output_path: Path) -> None:
        with self.handler(output_path, "w", **self.kwargs) as file:
            for distance in self.fixed:
                file.write(distance)
        assert_eq_files(output_path, self.fixed_path)


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
        x = Sequence("idx", self.seq_x)
        y = Sequence("idy", self.seq_y)
        r = self.metric.calculate(x, y)
        assert r.metric == self.metric
        assert r.x.id == "idx"
        assert r.y.id == "idy"
        if isinstance(r.d, float):
            assert abs(r.d - self.d) <= self.precision
        else:
            assert r.d == self.d
            assert r.d is None


class MetricFileTest(NamedTuple):
    file: str
    precision: float

    def get_metric_tests(self) -> Iterator[MetricTest]:
        path = TEST_DATA_DIR / self.file
        with DistanceHandler.Linear(path, "r") as file:
            for d in file:
                yield MetricTest(d.metric, d.x.id, d.y.id, d.d, self.precision)


def distances_simple() -> Distances:
    metric = DistanceMetric.Uncorrected()
    return Distances(
        [
            Distance(metric, Sequence("id1", None), Sequence("id2", None), 0.1),
            Distance(metric, Sequence("id1", None), Sequence("id3", None), 0.2),
            Distance(metric, Sequence("id1", None), Sequence("id4", None), 0.3),
        ]
    )


def distances_multiple() -> Distances:
    return Distances(
        [
            Distance(
                DistanceMetric.Uncorrected(),
                Sequence("id1", None),
                Sequence("id2", None),
                0.11,
            ),
            Distance(
                DistanceMetric.UncorrectedWithGaps(),
                Sequence("id1", None),
                Sequence("id2", None),
                0.12,
            ),
            Distance(
                DistanceMetric.JukesCantor(),
                Sequence("id1", None),
                Sequence("id2", None),
                0.13,
            ),
            Distance(
                DistanceMetric.Kimura2P(),
                Sequence("id1", None),
                Sequence("id2", None),
                0.14,
            ),
            Distance(
                DistanceMetric.NCD(), Sequence("id1", None), Sequence("id2", None), 0.15
            ),
            Distance(
                DistanceMetric.BBC(0),
                Sequence("id1", None),
                Sequence("id2", None),
                0.16,
            ),
            Distance(
                DistanceMetric.Uncorrected(),
                Sequence("id1", None),
                Sequence("id3", None),
                0.21,
            ),
            Distance(
                DistanceMetric.UncorrectedWithGaps(),
                Sequence("id1", None),
                Sequence("id3", None),
                0.22,
            ),
            Distance(
                DistanceMetric.JukesCantor(),
                Sequence("id1", None),
                Sequence("id3", None),
                0.23,
            ),
            Distance(
                DistanceMetric.Kimura2P(),
                Sequence("id1", None),
                Sequence("id3", None),
                0.24,
            ),
            Distance(
                DistanceMetric.NCD(), Sequence("id1", None), Sequence("id3", None), 0.25
            ),
            Distance(
                DistanceMetric.BBC(0),
                Sequence("id1", None),
                Sequence("id3", None),
                0.26,
            ),
            Distance(
                DistanceMetric.Uncorrected(),
                Sequence("id1", None),
                Sequence("id4", None),
                0.31,
            ),
            Distance(
                DistanceMetric.UncorrectedWithGaps(),
                Sequence("id1", None),
                Sequence("id4", None),
                0.32,
            ),
            Distance(
                DistanceMetric.JukesCantor(),
                Sequence("id1", None),
                Sequence("id4", None),
                0.33,
            ),
            Distance(
                DistanceMetric.Kimura2P(),
                Sequence("id1", None),
                Sequence("id4", None),
                0.34,
            ),
            Distance(
                DistanceMetric.NCD(), Sequence("id1", None), Sequence("id4", None), 0.35
            ),
            Distance(
                DistanceMetric.BBC(0),
                Sequence("id1", None),
                Sequence("id4", None),
                0.36,
            ),
        ]
    )


def distances_square() -> Distances:
    metric = DistanceMetric.Uncorrected()
    return Distances(
        [
            Distance(metric, Sequence("id1", None), Sequence("id1", None), 0.0),
            Distance(metric, Sequence("id1", None), Sequence("id2", None), 0.1),
            Distance(metric, Sequence("id1", None), Sequence("id3", None), 0.2),
            Distance(metric, Sequence("id2", None), Sequence("id1", None), 0.1),
            Distance(metric, Sequence("id2", None), Sequence("id2", None), 0.0),
            Distance(metric, Sequence("id2", None), Sequence("id3", None), 0.3),
            Distance(metric, Sequence("id3", None), Sequence("id1", None), 0.2),
            Distance(metric, Sequence("id3", None), Sequence("id2", None), 0.3),
            Distance(metric, Sequence("id3", None), Sequence("id3", None), 0.0),
        ]
    )


def distances_square_unknown() -> Distances:
    metric = DistanceMetric.Unknown()
    return Distances(
        [Distance(metric, dis.x, dis.y, dis.d) for dis in distances_square()]
    )


def distances_rectangle() -> Distances:
    metric = DistanceMetric.Uncorrected()
    return Distances(
        [
            Distance(metric, Sequence("id1", None), Sequence("id4", None), 0.14),
            Distance(metric, Sequence("id1", None), Sequence("id5", None), 0.15),
            Distance(metric, Sequence("id1", None), Sequence("id6", None), 0.16),
            Distance(metric, Sequence("id1", None), Sequence("id7", None), 0.17),
            Distance(metric, Sequence("id1", None), Sequence("id8", None), 0.18),
            Distance(metric, Sequence("id1", None), Sequence("id9", None), 0.19),
            Distance(metric, Sequence("id2", None), Sequence("id4", None), 0.24),
            Distance(metric, Sequence("id2", None), Sequence("id5", None), 0.25),
            Distance(metric, Sequence("id2", None), Sequence("id6", None), 0.26),
            Distance(metric, Sequence("id2", None), Sequence("id7", None), 0.27),
            Distance(metric, Sequence("id2", None), Sequence("id8", None), 0.28),
            Distance(metric, Sequence("id2", None), Sequence("id9", None), 0.29),
            Distance(metric, Sequence("id3", None), Sequence("id4", None), 0.34),
            Distance(metric, Sequence("id3", None), Sequence("id5", None), 0.35),
            Distance(metric, Sequence("id3", None), Sequence("id6", None), 0.36),
            Distance(metric, Sequence("id3", None), Sequence("id7", None), 0.37),
            Distance(metric, Sequence("id3", None), Sequence("id8", None), 0.38),
            Distance(metric, Sequence("id3", None), Sequence("id9", None), 0.39),
        ]
    )


def distances_missing() -> Distances:
    metric = DistanceMetric.Uncorrected()
    return Distances(
        [
            Distance(metric, Sequence("id1", None), Sequence("id1", None), 0.0),
            Distance(metric, Sequence("id1", None), Sequence("id2", None), None),
            Distance(metric, Sequence("id2", None), Sequence("id1", None), None),
            Distance(metric, Sequence("id2", None), Sequence("id2", None), 0.0),
        ]
    )


def distances_extras() -> Distances:
    return Distances(
        [
            Distance(
                DistanceMetric.Uncorrected(),
                Sequence("query1", None, dict(voucher="K")),
                Sequence("reference1", None, dict(voucher="X", organism="A")),
                0.11,
            ),
            Distance(
                DistanceMetric.UncorrectedWithGaps(),
                Sequence("query1", None, dict(voucher="K")),
                Sequence("reference1", None, dict(voucher="X", organism="A")),
                0.12,
            ),
            Distance(
                DistanceMetric.JukesCantor(),
                Sequence("query1", None, dict(voucher="K")),
                Sequence("reference1", None, dict(voucher="X", organism="A")),
                0.13,
            ),
            Distance(
                DistanceMetric.Kimura2P(),
                Sequence("query1", None, dict(voucher="K")),
                Sequence("reference1", None, dict(voucher="X", organism="A")),
                0.14,
            ),
            Distance(
                DistanceMetric.Uncorrected(),
                Sequence("query1", None, dict(voucher="K")),
                Sequence("reference2", None, dict(voucher="Y", organism="B")),
                0.21,
            ),
            Distance(
                DistanceMetric.UncorrectedWithGaps(),
                Sequence("query1", None, dict(voucher="K")),
                Sequence("reference2", None, dict(voucher="Y", organism="B")),
                0.22,
            ),
            Distance(
                DistanceMetric.JukesCantor(),
                Sequence("query1", None, dict(voucher="K")),
                Sequence("reference2", None, dict(voucher="Y", organism="B")),
                0.23,
            ),
            Distance(
                DistanceMetric.Kimura2P(),
                Sequence("query1", None, dict(voucher="K")),
                Sequence("reference2", None, dict(voucher="Y", organism="B")),
                0.24,
            ),
            Distance(
                DistanceMetric.Uncorrected(),
                Sequence("query2", None, dict(voucher="L")),
                Sequence("reference3", None, dict(voucher="Z", organism="C")),
                0.31,
            ),
            Distance(
                DistanceMetric.UncorrectedWithGaps(),
                Sequence("query2", None, dict(voucher="L")),
                Sequence("reference3", None, dict(voucher="Z", organism="C")),
                0.32,
            ),
            Distance(
                DistanceMetric.JukesCantor(),
                Sequence("query2", None, dict(voucher="L")),
                Sequence("reference3", None, dict(voucher="Z", organism="C")),
                0.33,
            ),
            Distance(
                DistanceMetric.Kimura2P(),
                Sequence("query2", None, dict(voucher="L")),
                Sequence("reference3", None, dict(voucher="Z", organism="C")),
                None,
            ),
        ]
    )


def distances_empty() -> Distances:
    return Distances([])


read_tests = [
    ReadTest(distances_simple, "simple.linear", DistanceHandler.Linear),
    ReadTest(distances_multiple, "multiple.linear", DistanceHandler.Linear),
    ReadTest(distances_missing, "missing.linear", DistanceHandler.Linear),
    ReadTest(distances_empty, "empty", DistanceHandler.Linear),
    ReadTest(distances_square_unknown, "square.matrix", DistanceHandler.Matrix),
    ReadTest(distances_empty, "empty", DistanceHandler.Matrix),
    ReadTest(
        distances_square,
        "square.matrix",
        DistanceHandler.Matrix,
        dict(metric=DistanceMetric.Uncorrected()),
    ),
    ReadTest(
        distances_rectangle,
        "rectangle.matrix",
        DistanceHandler.Matrix,
        dict(metric=DistanceMetric.Uncorrected()),
    ),
    ReadTest(
        distances_missing,
        "missing.matrix",
        DistanceHandler.Matrix,
        dict(metric=DistanceMetric.Uncorrected()),
    ),
    ReadTest(
        distances_extras,
        "extras.tsv",
        DistanceHandler.Linear.WithExtras,
        dict(idxHeader="seqid", idyHeader="id", tagX="_x", tagY="_y"),
    ),
    ReadTest(
        distances_extras,
        "extras.tsv",
        DistanceHandler.Linear.WithExtras,
        dict(idxColumn=0, idyColumn=2, tagX="_x", tagY="_y"),
    ),
    ReadTest(distances_empty, "empty", DistanceHandler.Linear.WithExtras),
]


write_tests = [
    WriteTest(
        distances_simple,
        "simple.linear",
        DistanceHandler.Linear,
        dict(formatter="{:.1f}"),
    ),
    WriteTest(
        distances_multiple,
        "multiple.linear",
        DistanceHandler.Linear,
        dict(formatter="{:.2f}"),
    ),
    WriteTest(
        distances_missing,
        "missing.linear",
        DistanceHandler.Linear,
        dict(formatter="{:.1f}"),
    ),
    WriteTest(
        distances_empty,
        "empty",
        DistanceHandler.Linear,
        dict(formatter="{:.1f}"),
    ),
    WriteTest(
        distances_square,
        "square.matrix",
        DistanceHandler.Matrix,
        dict(formatter="{:.1f}"),
    ),
    WriteTest(
        distances_rectangle,
        "rectangle.matrix",
        DistanceHandler.Matrix,
        dict(formatter="{:.2f}"),
    ),
    WriteTest(
        distances_missing,
        "missing.matrix",
        DistanceHandler.Matrix,
        dict(formatter="{:.1f}"),
    ),
    WriteTest(
        distances_empty,
        "empty",
        DistanceHandler.Matrix,
        dict(formatter="{:.1f}"),
    ),
    WriteTest(
        distances_missing,
        "missing.formatted.linear",
        DistanceHandler.Linear,
        dict(formatter="{:.2e}", missing="nan"),
    ),
    WriteTest(
        distances_missing,
        "missing.formatted.matrix",
        DistanceHandler.Matrix,
        dict(formatter="{:.2e}", missing="nan"),
    ),
    WriteTest(
        distances_extras,
        "extras.tsv",
        DistanceHandler.Linear.WithExtras,
        dict(
            idxHeader="seqid", idyHeader="id", tagX="_x", tagY="_y", formatter="{:.2f}"
        ),
    ),
    WriteTest(
        distances_missing,
        "missing.formatted.linear",
        DistanceHandler.Linear.WithExtras,
        dict(
            idxHeader="idx",
            idyHeader="idy",
            tagX="",
            tagY="",
            formatter="{:.2e}",
            missing="nan",
        ),
    ),
    WriteTest(
        distances_empty,
        "empty",
        DistanceHandler.Linear.WithExtras,
        dict(formatter="{:.1f}"),
    ),
]


label_tests = [
    LabelTest(DistanceMetric.Uncorrected(), "p"),
    LabelTest(DistanceMetric.UncorrectedWithGaps(), "p-gaps"),
    LabelTest(DistanceMetric.JukesCantor(), "jc"),
    LabelTest(DistanceMetric.Kimura2P(), "k2p"),
    LabelTest(DistanceMetric.NCD(), "ncd"),
    LabelTest(DistanceMetric.NCD(), "ncd"),
    LabelTest(DistanceMetric.BBC(0), "bbc(0)"),
    LabelTest(DistanceMetric.BBC(1), "bbc(1)"),
]


metric_tests = [
    MetricTest(DistanceMetric.Uncorrected(), "gg-ccnccta", "ggaccaccaa", 1.0 / 8.0),
    MetricTest(
        DistanceMetric.UncorrectedWithGaps(), "gg-ccnccta", "ggaccaccaa", 2.0 / 9.0
    ),
    MetricTest(DistanceMetric.Uncorrected(), "---", "nnn", None),
]


metric_file_tests = [
    MetricFileTest("metrics.tsv", 0.00051),
]


@pytest.mark.parametrize("test", read_tests)
def test_read_distances(test: ReadTest) -> None:
    test.validate()


@pytest.mark.parametrize("test", write_tests)
def test_write_distances(test: WriteTest, tmp_path: Path) -> None:
    output_path = test.get_output_path(tmp_path)
    test.validate(output_path)


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
        print(a.args[0], "\n", file=stderr)
    assert len(stack) == 0
