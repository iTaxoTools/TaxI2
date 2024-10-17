from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest
from utility import assert_eq_files

from itaxotools.taxi2.sequences import Sequence, SequenceHandler, Sequences

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
    fixture: Callable[[], Sequences]
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
        with self.handler(output_path, "w", **self.kwargs) as file:
            for sequence in self.fixed:
                file.write(sequence)
        assert_eq_files(output_path, self.fixed_path)


def sequences_simple() -> Sequences:
    return Sequences(
        [
            Sequence("id1", "ATC"),
            Sequence("id2", "ATG"),
            Sequence("id3", "ATA"),
        ]
    )


def sequences_organism() -> Sequences:
    return Sequences(
        [
            Sequence("id1", "ATC", {"organism": "X"}),
            Sequence("id2", "ATG", {"organism": "Y"}),
            Sequence("id3", "ATA", {"organism": "Z"}),
        ]
    )


def sequences_headers() -> Sequences:
    return Sequences(
        [
            Sequence("id1", "ATC", {"voucher": "X"}),
            Sequence("id2", "ATG", {"voucher": "Y"}),
            Sequence("id3", "ATA", {"voucher": "Z"}),
        ]
    )


def sequences_alleles() -> Sequences:
    return Sequences(
        [
            Sequence("id1", "ATC", {"allele": "a", "species": "X"}),
            Sequence("id1", "ATC", {"allele": "b", "species": "X"}),
            Sequence("id2", "ATG", {"allele": "a", "species": "Y"}),
            Sequence("id2", "ATG", {"allele": "b", "species": "Y"}),
            Sequence("id3", "ATA", {"allele": "a", "species": "Z"}),
            Sequence("id3", "ATA", {"allele": "b", "species": "Z"}),
        ]
    )


def sequences_quality() -> Sequences:
    return Sequences(
        [
            Sequence("id1", "ATC", {"quality": "!''"}),
            Sequence("id2", "ATG", {"quality": "(*)"}),
            Sequence("id3", "ATA", {"quality": "CF>"}),
        ]
    )


def sequences_empty() -> Sequences:
    return Sequences([])


@pytest.mark.parametrize(
    "test",
    [
        ReadTest(sequences_simple, "simple.fas", SequenceHandler.Fasta),
        ReadTest(sequences_simple, "simple.ali", SequenceHandler.Ali),
        ReadTest(sequences_simple, "simple.multi.fas", SequenceHandler.Fasta),
        ReadTest(sequences_simple, "simple.gbk", SequenceHandler.Genbank),
        ReadTest(sequences_simple, "simple.tsv", SequenceHandler.Tabfile),
        ReadTest(sequences_simple, "simple.xlsx", SequenceHandler.Excel),
        ReadTest(sequences_quality, "quality.fq", SequenceHandler.FastQ),
        ReadTest(
            sequences_headers,
            "headers.tsv",
            SequenceHandler.Tabfile,
            dict(idHeader="seqid", seqHeader="sequences"),
        ),
        ReadTest(
            sequences_headers,
            "headers.xlsx",
            SequenceHandler.Excel,
            dict(idHeader="seqid", seqHeader="sequences"),
        ),
        ReadTest(
            sequences_headers,
            "species.fas",
            SequenceHandler.Fasta,
            dict(parse_organism=True, organism_tag="voucher", organism_separator="|"),
        ),
        ReadTest(
            sequences_organism,
            "species.fas",
            SequenceHandler.Fasta,
            dict(parse_organism=True),
        ),
        ReadTest(
            sequences_organism,
            "species.dot.fas",
            SequenceHandler.Fasta,
            dict(parse_organism=True, organism_separator="."),
        ),
        ReadTest(sequences_empty, "empty", SequenceHandler.Fasta),
        ReadTest(sequences_empty, "empty", SequenceHandler.FastQ),
        ReadTest(sequences_empty, "empty", SequenceHandler.Ali),
        ReadTest(
            sequences_empty,
            "empty.tsv",
            SequenceHandler.Tabfile,
            dict(idHeader="seqid", seqHeader="sequences"),
        ),
    ],
)
def test_read_sequences(test: ReadTest) -> None:
    test.validate()


@pytest.mark.parametrize(
    "test",
    [
        WriteTest(sequences_simple, "simple.tsv", SequenceHandler.Tabfile),
        WriteTest(
            sequences_headers,
            "headers.tsv",
            SequenceHandler.Tabfile,
            dict(idHeader="seqid", seqHeader="sequences"),
        ),
        WriteTest(sequences_simple, "simple.fas", SequenceHandler.Fasta),
        WriteTest(
            sequences_simple,
            "simple.width.fas",
            SequenceHandler.Fasta,
            dict(line_width=2),
        ),
        WriteTest(
            sequences_organism,
            "species.fas",
            SequenceHandler.Fasta,
            dict(write_organism=True),
        ),
        WriteTest(
            sequences_organism,
            "species.dot.fas",
            SequenceHandler.Fasta,
            dict(write_organism=True, organism_separator="."),
        ),
        WriteTest(
            sequences_headers,
            "species.fas",
            SequenceHandler.Fasta,
            dict(write_organism=True, organism_tag="voucher", organism_separator="|"),
        ),
        WriteTest(
            sequences_alleles,
            "alleles.concat.fas",
            SequenceHandler.Fasta,
            dict(write_organism=False, concatenate_extras=["species", "allele"]),
        ),
        WriteTest(
            sequences_alleles,
            "alleles.plain.fas",
            SequenceHandler.Fasta,
            dict(write_organism=False, concatenate_extras=["allele"]),
        ),
        WriteTest(
            sequences_alleles,
            "alleles.species.fas",
            SequenceHandler.Fasta,
            dict(
                write_organism=True,
                organism_separator="|",
                organism_tag="species",
                concatenate_extras=["allele"],
            ),
        ),
    ],
)
def test_write_sequences(test: WriteTest, tmp_path: Path) -> None:
    output_path = test.get_output_path(tmp_path)
    test.validate(output_path)


def test_sequence_sanitized_id_with_extras():
    sequence = Sequence("samplé1", "ACGT", {"allele": "a"})
    id = sequence.get_sanitized_id_with_extras()
    assert id == "sample1_a"
