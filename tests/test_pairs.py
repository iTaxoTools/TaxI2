from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest
from utility import assert_eq_files

from itaxotools.taxi2.pairs import SequencePair, SequencePairHandler, SequencePairs
from itaxotools.taxi2.sequences import Sequence, Sequences

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class ReadTest(NamedTuple):
    fixture: Callable[[], SequencePairs]
    input: str
    handler: SequencePairHandler
    kwargs: dict = {}

    @property
    def input_path(self) -> Path:
        return TEST_DATA_DIR / self.input

    @property
    def fixed(self) -> SequencePairs:
        return self.fixture()

    def validate(self):
        pairs = SequencePairs.fromPath(self.input_path, self.handler, **self.kwargs)
        generated_list = list(pairs)
        fixed_list = list(self.fixed)
        assert len(fixed_list) == len(generated_list)
        for pair in fixed_list:
            assert pair in generated_list


class WriteTest(NamedTuple):
    fixture: Callable[[], SequencePairs]
    output: str
    handler: SequencePairHandler
    kwargs: dict = {}

    @property
    def fixed_path(self) -> Path:
        return TEST_DATA_DIR / self.output

    @property
    def fixed(self) -> SequencePairs:
        return self.fixture()

    def validate(self, tmp_path: Path) -> None:
        output_path = tmp_path / self.output
        print("!!", output_path)
        with self.handler(output_path, "w", **self.kwargs) as file:
            for pair in self.fixed:
                file.write(pair)
        assert_eq_files(output_path, self.fixed_path, ignore=r"\n")


def pairs_simple() -> SequencePairs:
    return SequencePairs(
        [
            SequencePair(
                Sequence("id1", "ATC-"),
                Sequence("id2", "ATG-"),
            ),
            SequencePair(
                Sequence("id1", "ATC-"),
                Sequence("id3", "-TAA"),
            ),
            SequencePair(
                Sequence("id2", "ATG-"),
                Sequence("id3", "-TAA"),
            ),
        ]
    )


read_tests = [
    ReadTest(pairs_simple, "simple.tsv", SequencePairHandler.Tabfile),
    ReadTest(pairs_simple, "simple.formatted", SequencePairHandler.Formatted),
]


write_tests = [
    WriteTest(pairs_simple, "simple.tsv", SequencePairHandler.Tabfile),
    WriteTest(pairs_simple, "simple.formatted", SequencePairHandler.Formatted),
]


def test_pairs_from_product() -> None:
    xs = [
        Sequence("id1", "ATC"),
        Sequence("id2", "ATG"),
    ]
    ys = [
        Sequence("id3", "TAA"),
        Sequence("id4", "TAC"),
        Sequence("id5", "TAG"),
    ]
    ts = [
        SequencePair(xs[0], ys[0]),
        SequencePair(xs[0], ys[1]),
        SequencePair(xs[0], ys[2]),
        SequencePair(xs[1], ys[0]),
        SequencePair(xs[1], ys[1]),
        SequencePair(xs[1], ys[2]),
    ]
    ps = SequencePairs.fromProduct(Sequences(xs), Sequences(ys))
    for p, t in zip(ps, ts):
        assert p == t


@pytest.mark.parametrize("test", read_tests)
def test_read_pairs(test: ReadTest) -> None:
    test.validate()


@pytest.mark.parametrize("test", write_tests)
def test_write_pairs(test: WriteTest, tmp_path: Path) -> None:
    test.validate(tmp_path)
