from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest
from utility import assert_eq_files

from itaxotools.taxi3.sequences import Sequence, Sequences, SequenceHandler

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class ReadTest(NamedTuple):
    fixture: Callable[[], Sequences]
    input: str
    handler: SequenceHandler
    kwargs: dict = {}

    @property
    def input_path(self) -> Path:
        return TEST_DATA_DIR / self.input

    @property
    def fixed(self) -> Sequences:
        return self.fixture()

    def validate(self) -> None:
        sequences = Sequences.fromPath(self.input_path, self.handler, **self.kwargs)
        generated_list = list(sequences)
        fixed_list = list(self.fixed)
        assert len(fixed_list) == len(generated_list)
        for sequence in fixed_list:
            assert sequence in generated_list


class WriteTest(NamedTuple):
    fixture: Callable[[], Items]
    output: str
    handler: SequenceHandler
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
        with self.handler(output_path, 'w', **self.kwargs) as file:
            for sequence in self.fixed:
                file.write(sequence)
        assert_eq_files(output_path, self.fixed_path)


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


def sequences_empty() -> Sequences:
    return Sequences([])


@pytest.mark.parametrize(
    "test", [
    ReadTest(sequences_simple, 'simple.fas', SequenceHandler.Fasta),
    ReadTest(sequences_simple, 'simple.multi.fas', SequenceHandler.Fasta),
    ReadTest(sequences_simple, 'simple.gbk', SequenceHandler.Genbank),
    ReadTest(sequences_simple, 'simple.tsv', SequenceHandler.Tabfile),
    ReadTest(sequences_simple, 'simple.xlsx', SequenceHandler.Excel),
    ReadTest(sequences_headers, 'headers.tsv', SequenceHandler.Tabfile, dict(idHeader='seqid', seqHeader='sequences')),
    ReadTest(sequences_headers, 'headers.xlsx', SequenceHandler.Excel, dict(idHeader='seqid', seqHeader='sequences')),
    ReadTest(sequences_empty, 'empty.tsv', SequenceHandler.Tabfile, dict(idHeader='seqid', seqHeader='sequences')),
])
def test_read_sequences(test: ReadTest) -> None:
    test.validate()


@pytest.mark.parametrize(
    "test", [
    WriteTest(sequences_simple, 'simple.tsv', SequenceHandler.Tabfile),
    WriteTest(sequences_headers, 'headers.tsv', SequenceHandler.Tabfile, dict(idHeader='seqid', seqHeader='sequences')),
])
def test_write_sequences(test: WriteTest, tmp_path: Path) -> None:
    output_path = test.get_output_path(tmp_path)
    test.validate(output_path)
