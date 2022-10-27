from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest
from utility import assert_eq_files

from itaxotools.taxi3.pairs import (
    SequencePair, SequencePairFile, SequencePairs)
from itaxotools.taxi3.sequences import Sequence

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class ReadTest(NamedTuple):
    fixture: Callable[[], SequencePairs]
    input: str
    file: SequencePairFile
    kwargs: dict = {}

    def validate(self, generated: SequencePairs):
        fixture_list = list(self.fixture())
        generated_list = list(generated)
        for pair in fixture_list:
            assert pair in generated_list


class WriteTest(NamedTuple):
    fixture: Callable[[], SequencePairs]
    output: str
    file: SequencePairFile

    def generate(self) -> SequencePairs:
        return self.fixture()


def pairs_simple() -> SequencePairs:
    return SequencePairs([
        SequencePair(
            Sequence('id1', 'ATC-'),
            Sequence('id2', 'ATG-'),
        ),
        SequencePair(
            Sequence('id1', 'ATC-'),
            Sequence('id3', '-TAA'),
        ),
        SequencePair(
            Sequence('id2', 'ATG-'),
            Sequence('id3', '-TAA'),
        ),
    ])


read_tests = [
    ReadTest(pairs_simple, 'simple.tsv', SequencePairFile.Tabfile),
    ReadTest(pairs_simple, 'simple.formatted', SequencePairFile.Formatted),
]


write_tests = [
    WriteTest(pairs_simple, 'simple.tsv', SequencePairFile.Tabfile),
    WriteTest(pairs_simple, 'simple.formatted', SequencePairFile.Formatted),
]


@pytest.mark.xfail
@pytest.mark.parametrize("test", read_tests)
def test_read_pairs(test: ReadTest) -> None:
    input_path = TEST_DATA_DIR / test.input
    pairs = test.file(input_path).read(**test.kwargs)
    test.validate(pairs)


@pytest.mark.xfail
@pytest.mark.parametrize("test", write_tests)
def test_write_pairs(test: WriteTest, tmp_path: Path) -> None:
    fixed_path = TEST_DATA_DIR / test.output
    output_path = tmp_path / test.output
    pairs = iter(test.generate())
    test.file(output_path).write(pairs)
    assert_eq_files(output_path, fixed_path, ignore=r'\n')
