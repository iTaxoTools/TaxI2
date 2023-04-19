from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple
from os.path import getsize

import pytest
from utility import assert_eq_files

from itaxotools.taxi2.files import FileFormat, FileInfo

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
        assert result.size == getsize(str(self.input_path))
        for info in self.infos:
            assert getattr(result, info) == self.infos[info]


@pytest.mark.parametrize(
    "test", [
    IdentifyTest(FileFormat.Fasta, 'simple.fasta'),
    IdentifyTest(FileFormat.Tabfile, 'simple.tsv'),
    IdentifyTest(FileFormat.Spart, 'simple.spart'),
    IdentifyTest(FileFormat.Spart, 'simple.xml'),
    IdentifyTest(FileFormat.Unknown, 'empty.txt'),
])
def test_identify_file(test: IdentifyTest) -> None:
    test.validate()


@pytest.mark.parametrize(
    "test", [
    InfoTest('simple.fasta', dict(format=FileFormat.Fasta, has_subsets=False)),
    InfoTest('species.fasta', dict(format=FileFormat.Fasta, has_subsets=True)),
    InfoTest('simple.fasta', dict(format=FileFormat.Fasta, has_subsets=False)),
    InfoTest('simple.tsv', dict(format=FileFormat.Tabfile, headers=['id', 'seq'])),
    InfoTest('full.tsv', dict(
        format=FileFormat.Tabfile,
        headers=['seqid', 'voucher', 'organism', 'genus', 'species', 'sequence'],
        header_individuals='seqid',
        header_sequences='sequence',
        header_organism='organism',
        header_species='species',
        header_genus='genus',
        )),
    InfoTest('binomen.tsv', dict(
        format=FileFormat.Tabfile,
        headers=['seqid', 'voucher', 'species', 'sequence'],
        header_individuals='seqid',
        header_sequences='sequence',
        header_organism='species',
        header_species=None,
        header_genus=None,
        )),
    InfoTest('simple.spart', dict(format=FileFormat.Spart, spartitions=['spartition_1', 'spartition_2'], is_matricial=True, is_xml=False)),
    InfoTest('simple.xml', dict(format=FileFormat.Spart, spartitions=['spartition_1'], is_matricial=False, is_xml=True)),
    InfoTest('empty.txt', dict(format=FileFormat.Unknown)),
])
def test_get_file_info(test: InfoTest) -> None:
    test.validate()
