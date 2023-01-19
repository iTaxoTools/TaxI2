from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest
from utility import assert_eq_files

from itaxotools.taxi3.files import FileFormat, FileInfo

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class IdentifyTest(NamedTuple):
    format: FileFormat
    input: str

    @property
    def input_path(self) -> Path:
        return TEST_DATA_DIR / self.input

    def validate(self) -> None:
        result = FileFormat.identify(self.input_path)
        assert result == self.format


class InfoTest(NamedTuple):
    input: str
    infos: dict[str, object]

    @property
    def input_path(self) -> Path:
        return TEST_DATA_DIR / self.input

    def validate(self) -> None:
        result = FileInfo.from_path(self.input_path)
        assert result.path == self.input_path
        for info in self.infos:
            assert getattr(result, info) == self.infos[info]


@pytest.mark.parametrize(
    "test", [
    IdentifyTest(FileFormat.Fasta, 'simple.fasta'),
    IdentifyTest(FileFormat.Tabfile, 'simple.tsv'),
    IdentifyTest(FileFormat.Unknown, 'empty.txt'),
])
def test_identify_file(test: IdentifyTest) -> None:
    test.validate()


@pytest.mark.parametrize(
    "test", [
    InfoTest('simple.fasta', dict(format=FileFormat.Fasta, size=15)),
    InfoTest('simple.tsv', dict(format=FileFormat.Tabfile, size=21, headers=['id', 'seq'])),
    InfoTest('full.tsv', dict(
        format=FileFormat.Tabfile,
        headers=['seqid', 'voucher', 'organism', 'genus', 'species', 'sequence'],
        header_individuals='seqid',
        header_sequences='sequence',
        header_organism='organism',
        header_species='species',
        header_genus='genus',
        )),
    InfoTest('empty.txt', dict(format=FileFormat.Unknown, size=0)),
])
def test_get_file_info(test: InfoTest) -> None:
    test.validate()
