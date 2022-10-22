from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest

from itaxotools.taxi3.sequences import Sequence, Sequences, SequenceReader

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class ReadTest(NamedTuple):
    validator: Callable
    input: str
    reader: SequenceReader
    kwargs: dict


def validate_simple(sequences: Sequences):
    seq_list = list(sequences)
    assert len(seq_list) == 3
    assert Sequence('id1', 'ATC') in seq_list
    assert Sequence('id2', 'ATG') in seq_list
    assert Sequence('id3', 'ATA') in seq_list


read_tests = [
    ReadTest(validate_simple, 'simple.fas', SequenceReader.FastaReader, {}),
    ReadTest(validate_simple, 'simple.multi.fas', SequenceReader.FastaReader, {}),
    ReadTest(validate_simple, 'simple.gbk', SequenceReader.GenbankReader, {}),
    ReadTest(validate_simple, 'simple.tsv', SequenceReader.TabfileReader, {}),
    ReadTest(validate_simple, 'simple.headers.tsv', SequenceReader.TabfileReader, dict(idHeader='seqid', seqHeader='sequences')),
    ReadTest(validate_simple, 'simple.xlsx', SequenceReader.ExcelReader, {}),
    ReadTest(validate_simple, 'simple.headers.xlsx', SequenceReader.ExcelReader, dict(id='seqid', seq='sequence')),
]


@pytest.mark.parametrize("test", read_tests)
def test_read_sequences(test: ReadTest) -> None:
    input_path = TEST_DATA_DIR / test.input
    sequences = test.reader.read(input_path, **test.kwargs)
    test.validator(sequences)
