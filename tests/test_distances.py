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


def validate_simple(distances: Distances):
    dis_list = list(distances)
    assert len(dis_list) == 3
    assert Distance(DistanceMetric.Uncorrected(), 'id1', 'id2', 0.1) in dis_list
    assert Distance(DistanceMetric.Uncorrected(), 'id1', 'id3', 0.2) in dis_list
    assert Distance(DistanceMetric.Uncorrected(), 'id1', 'id4', 0.3) in dis_list


def generate_simple() -> Distances:
    distances = [
        Distance(DistanceMetric.Uncorrected(), 'id1', 'id2', 0.1),
        Distance(DistanceMetric.Uncorrected(), 'id1', 'id3', 0.2),
        Distance(DistanceMetric.Uncorrected(), 'id1', 'id4', 0.3),
    ]
    return Distances(lambda: iter(distances))


read_tests = [
    ReadTest(validate_simple, 'simple.tsv', DistanceFile.Linear),
]


write_tests = [
    WriteTest(generate_simple, 'simple.tsv', DistanceFile.Linear),
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
