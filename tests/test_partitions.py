from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest

from itaxotools.taxi2.pairs import SequencePairs
from itaxotools.taxi2.partitions import Classification, Partition, PartitionHandler

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class ReadTest(NamedTuple):
    fixture: Callable[[], SequencePairs]
    input: str
    handler: PartitionHandler
    kwargs: dict = {}

    @property
    def input_path(self) -> Path:
        return TEST_DATA_DIR / self.input

    @property
    def fixed(self) -> SequencePairs:
        return self.fixture()

    def validate(self):
        generated = Partition.fromPath(self.input_path, self.handler, **self.kwargs)
        assert generated == self.fixed


def spartition_simple() -> Partition:
    return {
        "sample1": "speciesA",
        "sample2": "speciesA",
        "sample3": "speciesA",
        "sample4": "speciesA",
        "sample5": "speciesB",
        "sample6": "speciesB",
        "sample7": "speciesC",
    }


def spartition_missing() -> Partition:
    return {
        "sample3": "speciesA",
        "sample4": "speciesA",
        "sample6": "speciesB",
        "sample7": "speciesC",
    }


def spartition_matricial() -> Partition:
    return {
        "sample1": "1",
        "sample2": "1",
        "sample3": "1",
        "sample4": "1",
        "sample5": "2",
        "sample6": "2",
        "sample7": "3",
    }


def spartition_genera() -> Partition:
    return {
        "sample1": "genusX",
        "sample2": "genusX",
        "sample3": "genusX",
        "sample4": "genusX",
        "sample5": "genusY",
        "sample6": "genusY",
        "sample7": "genusY",
    }


read_tests = [
    ReadTest(spartition_simple, "simple.tsv", PartitionHandler.Tabfile),
    ReadTest(
        spartition_simple,
        "extras.tsv",
        PartitionHandler.Tabfile,
        dict(idHeader="seqid", subHeader="organism"),
    ),
    ReadTest(
        spartition_genera,
        "genera.tsv",
        PartitionHandler.Tabfile,
        dict(
            filter=PartitionHandler.subset_first_word,
            idHeader="seqid",
            subHeader="organism",
        ),
    ),
    ReadTest(spartition_simple, "simple.xml", PartitionHandler.Spart),
    ReadTest(spartition_matricial, "simple.spart", PartitionHandler.Spart),
    ReadTest(spartition_simple, "simple.fas", PartitionHandler.Fasta),
    ReadTest(
        spartition_simple, "simple.dot.fas", PartitionHandler.Fasta, dict(separator=".")
    ),
    ReadTest(spartition_missing, "missing.fas", PartitionHandler.Fasta),
    ReadTest(
        spartition_genera,
        "genera.fas",
        PartitionHandler.Fasta,
        dict(filter=PartitionHandler.subset_first_word),
    ),
    ReadTest(
        spartition_simple,
        "genera.fas",
        PartitionHandler.Fasta,
        dict(filter=lambda x: Classification(x.individual, x.subset.split(" ")[1])),
    ),
]


@pytest.mark.parametrize("test", read_tests)
def test_read_pairs(test: ReadTest) -> None:
    test.validate()
