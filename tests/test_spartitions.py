from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest
from utility import assert_eq_files

from itaxotools.taxi3.spartitions import Spartition, SpartitionFile

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class ReadTest(NamedTuple):
    fixture: Callable[[], SequencePairs]
    input: str
    file: SpartitionFile
    kwargs: dict = {}

    def validate(self, generated: Spartition):
        fixture = self.fixture()
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


read_tests = [
    ReadTest(spartition_simple, 'simple.tsv', SpartitionFile.Tabfile, dict(idHeader='seqid', subsetHeader='organism')),
]


@pytest.mark.parametrize("test", read_tests)
def test_read_pairs(test: ReadTest) -> None:
    input_path = TEST_DATA_DIR / test.input
    spartition = test.file(input_path).read(**test.kwargs)
    test.validate(spartition)
