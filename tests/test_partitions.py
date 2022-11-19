from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest
from utility import assert_eq_files

from itaxotools.taxi3.partitions import Partition, PartitionFile

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class ReadTest(NamedTuple):
    fixture: Callable[[], SequencePairs]
    input: str
    file: PartitionFile
    kwargs: dict = {}

    def validate(self, generated: Partition):
        fixture = self.fixture()
        print(generated)
        assert fixture == generated


def spartition_simple() -> SequencePairs:
    return {
        'sample1': 'speciesA',
        'sample2': 'speciesA',
        'sample3': 'speciesA',
        'sample4': 'speciesA',
        'sample5': 'speciesB',
        'sample6': 'speciesB',
        'sample7': 'speciesC',
    }


def spartition_matricial() -> SequencePairs:
    return {
        'sample1': '1',
        'sample2': '1',
        'sample3': '1',
        'sample4': '1',
        'sample5': '2',
        'sample6': '2',
        'sample7': '3',
    }


read_tests = [
    ReadTest(spartition_simple, 'simple.tsv', PartitionFile.Tabfile, dict(idHeader='seqid', subsetHeader='organism')),
    ReadTest(spartition_simple, 'simple.xml', PartitionFile.Spart),
    ReadTest(spartition_matricial, 'simple.spart', PartitionFile.Spart),
]


@pytest.mark.parametrize("test", read_tests)
def test_read_pairs(test: ReadTest) -> None:
    input_path = TEST_DATA_DIR / test.input
    spartition = test.file(input_path).get(**test.kwargs)
    test.validate(spartition)
