from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest

from utility import assert_eq_files

from itaxotools.taxi3.distances import Distance, Distances, DistanceFile, DistanceMetric

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class ReadTest(NamedTuple):
    fixture: Callable[None, Distances]
    input: str
    file: DistanceFile
    kwargs: dict = {}

    def validate(self, generated: Distances):
        fixture_list = list(self.fixture())
        generated_list = list(generated)
        assert len(fixture_list) == len(generated_list)
        for distance in fixture_list:
            assert distance in generated_list


class WriteTest(NamedTuple):
    fixture: Callable[None, Distances]
    output: str
    file: DistanceFile

    def generate(self) -> Distances:
        return self.fixture()


def distances_simple() -> Distances:
    metric = DistanceMetric.Uncorrected()
    return Distances([
        Distance(metric, 'id1', 'id2', 0.1),
        Distance(metric, 'id1', 'id3', 0.2),
        Distance(metric, 'id1', 'id4', 0.3),
    ])


def distances_multiple() -> Distances:
    return Distances([
        Distance(DistanceMetric.Uncorrected(), 'id1', 'id2', 0.11),
        Distance(DistanceMetric.UncorrectedWithGaps(), 'id1', 'id2', 0.12),
        Distance(DistanceMetric.JukesCantor(), 'id1', 'id2', 0.13),
        Distance(DistanceMetric.Kimura2P(), 'id1', 'id2', 0.14),
        Distance(DistanceMetric.NCD(0), 'id1', 'id2', 0.15),
        Distance(DistanceMetric.BBC(0), 'id1', 'id2', 0.16),

        Distance(DistanceMetric.Uncorrected(), 'id1', 'id3', 0.21),
        Distance(DistanceMetric.UncorrectedWithGaps(), 'id1', 'id3', 0.22),
        Distance(DistanceMetric.JukesCantor(), 'id1', 'id3', 0.23),
        Distance(DistanceMetric.Kimura2P(), 'id1', 'id3', 0.24),
        Distance(DistanceMetric.NCD(0), 'id1', 'id3', 0.25),
        Distance(DistanceMetric.BBC(0), 'id1', 'id3', 0.26),

        Distance(DistanceMetric.Uncorrected(), 'id1', 'id4', 0.31),
        Distance(DistanceMetric.UncorrectedWithGaps(), 'id1', 'id4', 0.32),
        Distance(DistanceMetric.JukesCantor(), 'id1', 'id4', 0.33),
        Distance(DistanceMetric.Kimura2P(), 'id1', 'id4', 0.34),
        Distance(DistanceMetric.NCD(0), 'id1', 'id4', 0.35),
        Distance(DistanceMetric.BBC(0), 'id1', 'id4', 0.36),
    ])


def distances_square() -> Distances:
    metric = DistanceMetric.Uncorrected()
    return Distances([
        Distance(metric, 'id1', 'id1', 0.0),
        Distance(metric, 'id1', 'id2', 0.1),
        Distance(metric, 'id1', 'id3', 0.2),

        Distance(metric, 'id2', 'id1', 0.1),
        Distance(metric, 'id2', 'id2', 0.0),
        Distance(metric, 'id2', 'id3', 0.3),

        Distance(metric, 'id3', 'id1', 0.2),
        Distance(metric, 'id3', 'id2', 0.3),
        Distance(metric, 'id3', 'id3', 0.0),
    ])


def distances_square_unknown() -> Distances:
    metric = DistanceMetric.Unknown()
    return Distances([
        Distance(metric, dis.idx, dis.idy, dis.d) for dis in distances_square()
    ])


def distances_rectangle() -> Distances:
    metric = DistanceMetric.Uncorrected()
    return Distances([
        Distance(metric, 'id1', 'id4', 0.14),
        Distance(metric, 'id1', 'id5', 0.15),
        Distance(metric, 'id1', 'id6', 0.16),
        Distance(metric, 'id1', 'id7', 0.17),
        Distance(metric, 'id1', 'id8', 0.18),
        Distance(metric, 'id1', 'id9', 0.19),

        Distance(metric, 'id2', 'id4', 0.24),
        Distance(metric, 'id2', 'id5', 0.25),
        Distance(metric, 'id2', 'id6', 0.26),
        Distance(metric, 'id2', 'id7', 0.27),
        Distance(metric, 'id2', 'id8', 0.28),
        Distance(metric, 'id2', 'id9', 0.29),

        Distance(metric, 'id3', 'id4', 0.34),
        Distance(metric, 'id3', 'id5', 0.35),
        Distance(metric, 'id3', 'id6', 0.36),
        Distance(metric, 'id3', 'id7', 0.37),
        Distance(metric, 'id3', 'id8', 0.38),
        Distance(metric, 'id3', 'id9', 0.39),
    ])


read_tests = [
    ReadTest(distances_simple, 'simple.tsv', DistanceFile.Linear),
    ReadTest(distances_multiple, 'multiple.tsv', DistanceFile.Linear),
    ReadTest(distances_square_unknown, 'square.tsv', DistanceFile.Matrix),
    ReadTest(distances_square, 'square.tsv', DistanceFile.Matrix, dict(metric=DistanceMetric.Uncorrected())),
    ReadTest(distances_rectangle, 'rectangle.tsv', DistanceFile.Matrix, dict(metric=DistanceMetric.Uncorrected())),
]


write_tests = [
    WriteTest(distances_simple, 'simple.tsv', DistanceFile.Linear),
    WriteTest(distances_multiple, 'multiple.tsv', DistanceFile.Linear),
    WriteTest(distances_square, 'square.tsv', DistanceFile.Matrix),
    WriteTest(distances_rectangle, 'rectangle.tsv', DistanceFile.Matrix),
]


@pytest.mark.parametrize("test", read_tests)
def test_read_distances(test: ReadTest) -> None:
    input_path = TEST_DATA_DIR / test.input
    distances = test.file(input_path).read(**test.kwargs)
    test.validate(distances)


@pytest.mark.parametrize("test", write_tests)
def test_write(test: WriteTest, tmp_path: Path) -> None:
    fixed_path = TEST_DATA_DIR / test.output
    output_path = tmp_path / test.output
    distances = test.generate()
    test.file(output_path).write(distances)
    assert_eq_files(output_path, fixed_path)
