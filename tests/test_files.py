from __future__ import annotations

from os.path import getsize
from pathlib import Path
from typing import NamedTuple

import pytest

from itaxotools.taxi2.file_types import FileFormat
from itaxotools.taxi2.files import get_info, identify_format

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class IdentifyTest(NamedTuple):
    format: FileFormat
    input: str

    @property
    def input_path(self) -> Path:
        return TEST_DATA_DIR / self.input

    def validate(self) -> None:
        result = identify_format(self.input_path)
        assert result == self.format


class InfoTest(NamedTuple):
    input: str
    infos: dict[str, object]

    @property
    def input_path(self) -> Path:
        return TEST_DATA_DIR / self.input

    def validate(self) -> None:
        result = get_info(self.input_path)
        assert result.path == self.input_path
        assert result.size == getsize(str(self.input_path))
        for info in self.infos:
            assert getattr(result, info) == self.infos[info]


@pytest.mark.parametrize(
    "test",
    [
        IdentifyTest(FileFormat.Fasta, "simple.fasta"),
        IdentifyTest(FileFormat.Tabfile, "simple.tsv"),
        IdentifyTest(FileFormat.Spart, "simple.spart"),
        IdentifyTest(FileFormat.Spart, "simple.xml"),
        IdentifyTest(FileFormat.Unknown, "empty.txt"),
    ],
)
def test_identify_file(test: IdentifyTest) -> None:
    test.validate()


@pytest.mark.parametrize(
    "test",
    [
        InfoTest("simple.fasta", dict(format=FileFormat.Fasta, has_subsets=False)),
        InfoTest("spaces.fasta", dict(format=FileFormat.Fasta, has_subsets=False)),
        InfoTest(
            "species.fasta",
            dict(format=FileFormat.Fasta, has_subsets=True, subset_separator="|"),
        ),
        InfoTest(
            "species.dot.fasta",
            dict(format=FileFormat.Fasta, has_subsets=True, subset_separator="."),
        ),
        InfoTest("simple.fasta", dict(format=FileFormat.Fasta, has_subsets=False)),
        InfoTest("simple.tsv", dict(format=FileFormat.Tabfile, headers=["id", "seq"])),
        InfoTest(
            "full.tsv",
            dict(
                format=FileFormat.Tabfile,
                headers=[
                    "seqid",
                    "voucher",
                    "organism",
                    "genus",
                    "species",
                    "sequence",
                ],
                header_individuals="seqid",
                header_sequences="sequence",
                header_organism="organism",
                header_species="species",
                header_genus="genus",
            ),
        ),
        InfoTest(
            "binomen.tsv",
            dict(
                format=FileFormat.Tabfile,
                headers=["seqid", "voucher", "species", "sequence"],
                header_individuals="seqid",
                header_sequences="sequence",
                header_organism="species",
                header_species=None,
                header_genus=None,
            ),
        ),
        InfoTest(
            "simple.spart",
            dict(
                format=FileFormat.Spart,
                spartitions=["spartition_1", "spartition_2"],
                is_matricial=True,
                is_xml=False,
            ),
        ),
        InfoTest(
            "simple.xml",
            dict(
                format=FileFormat.Spart,
                spartitions=["spartition_1"],
                is_matricial=False,
                is_xml=True,
            ),
        ),
        InfoTest(
            "simple.tree",
            dict(format=FileFormat.Newick, count=1, names={"A", "B", "C", "D"}),
        ),
        InfoTest(
            "complex.tree",
            dict(format=FileFormat.Newick, count=2, names={"A", "B", "C", "D", "E"}),
        ),
        InfoTest("empty.txt", dict(format=FileFormat.Unknown)),
    ],
)
def test_get_file_info(test: InfoTest) -> None:
    test.validate()
