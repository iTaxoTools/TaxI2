from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest

from itaxotools.taxi3.sequences import Sequence, SequenceFile, Sequences

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class ReadTest(NamedTuple):
    fixture: Callable[[], Sequences]
    input: str
    file: SequenceFile
    kwargs: dict = {}

    def validate(self, generated: Sequences):
        fixture_list = list(self.fixture())
        generated_list = list(generated)
        assert len(fixture_list) == len(generated_list)
        for sequence in fixture_list:
            assert sequence in generated_list


def sequences_simple() -> Sequences:
    return Sequences([
        Sequence('id1', 'ATC'),
        Sequence('id2', 'ATG'),
        Sequence('id3', 'ATA'),
    ])


def sequences_headers() -> Sequences:
    return Sequences([
        Sequence('id1', 'ATC', {'voucher': 'X'}),
        Sequence('id2', 'ATG', {'voucher': 'Y'}),
        Sequence('id3', 'ATA', {'voucher': 'Z'}),
    ])


read_tests = [
    ReadTest(sequences_simple, 'simple.fas', SequenceFile.Fasta),
    ReadTest(sequences_simple, 'simple.multi.fas', SequenceFile.Fasta),
    ReadTest(sequences_simple, 'simple.gbk', SequenceFile.Genbank),
    ReadTest(sequences_simple, 'simple.tsv', SequenceFile.Tabfile),
    ReadTest(sequences_simple, 'simple.xlsx', SequenceFile.Excel),
    ReadTest(sequences_headers, 'headers.tsv', SequenceFile.Tabfile, dict(idHeader='seqid', seqHeader='sequences')),
    ReadTest(sequences_headers, 'headers.xlsx', SequenceFile.Excel, dict(idHeader='seqid', seqHeader='sequences')),
]


@pytest.mark.parametrize("test", read_tests)
def test_read_sequences(test: ReadTest) -> None:
    input_path = TEST_DATA_DIR / test.input
    sequences = test.file(input_path).read(**test.kwargs)
    test.validate(sequences)
