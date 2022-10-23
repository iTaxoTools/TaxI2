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
    return Distances([
        Distance(DistanceMetric.Uncorrected(), 'id1', 'id2', 0.1),
        Distance(DistanceMetric.Uncorrected(), 'id1', 'id3', 0.2),
        Distance(DistanceMetric.Uncorrected(), 'id1', 'id4', 0.3),
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


read_tests = [
    ReadTest(distances_simple, 'simple.tsv', DistanceFile.Linear),
    ReadTest(distances_multiple, 'multiple.tsv', DistanceFile.Linear),
    ReadTest(distances_multiple, 'multiple.tsv', DistanceFile.Matrix),
]


write_tests = [
    WriteTest(distances_simple, 'simple.tsv', DistanceFile.Linear),
    WriteTest(distances_multiple, 'multiple.tsv', DistanceFile.Linear),
]


@pytest.mark.parametrize("test", read_tests)
def test_read_distances(test: ReadTest) -> None:
    input_path = TEST_DATA_DIR / test.input
    distances = test.file(input_path).read()
    test.validate(distances)


@pytest.mark.parametrize("test", write_tests)
def test_write(test: WriteTest, tmp_path: Path) -> None:
    fixed_path = TEST_DATA_DIR / test.output
    output_path = tmp_path / test.output
    distances = test.generate()
    test.file(output_path).write(distances)
    assert_eq_files(output_path, fixed_path)
