from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest

from utility import assert_eq_files

from itaxotools.taxi3.distances import Distance, Distances, DistanceFile, DistanceMetric

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class ReadTest(NamedTuple):
    validator: Callable
    input: str
    file: DistanceFile


class WriteTest(NamedTuple):
    generator: Callable
    output: str
    file: DistanceFile


def distances_multiple():
    return [
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
    ]


def distances_simple():
    return [
        Distance(DistanceMetric.Uncorrected(), 'id1', 'id2', 0.1),
        Distance(DistanceMetric.Uncorrected(), 'id1', 'id3', 0.2),
        Distance(DistanceMetric.Uncorrected(), 'id1', 'id4', 0.3),
    ]


def validate_simple(distances: Distances):
    test_list = distances_simple()
    dis_list = list(distances)
    assert len(dis_list) == len(test_list)
    for distance in test_list:
        assert distance in dis_list


def validate_multiple(distances: Distances):
    test_list = distances_multiple()
    dis_list = list(distances)
    assert len(dis_list) == len(test_list)

    for distance in test_list:
        assert distance in dis_list


def generate_simple() -> Distances:
    distances = distances_simple()
    return Distances(lambda: iter(distances))


def generate_multiple() -> Distances:
    distances = distances_multiple()
    return Distances(lambda: iter(distances))


read_tests = [
    ReadTest(validate_simple, 'simple.tsv', DistanceFile.Linear),
    ReadTest(validate_simple, 'multiple.tsv', DistanceFile.Linear),
]


write_tests = [
    WriteTest(generate_simple, 'simple.tsv', DistanceFile.Linear),
    WriteTest(generate_multiple, 'multiple.tsv', DistanceFile.Linear),
]


@pytest.mark.parametrize("test", read_tests)
def test_read_distances(test: ReadTest) -> None:
    input_path = TEST_DATA_DIR / test.input
    distances = test.file(input_path).read()
    test.validator(distances)


@pytest.mark.parametrize("test", write_tests)
def test_write(test: WriteTest, tmp_path: Path) -> None:
    output_path = TEST_DATA_DIR / test.output
    test_path = tmp_path / test.output
    distances = test.generator()
    test.file(test_path).write(distances)
    assert_eq_files(test_path, output_path, test.case_sensitive)
