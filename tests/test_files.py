from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest
from utility import assert_eq_files

from itaxotools.taxi3.files import FileFormat

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


@pytest.mark.parametrize(
    "test", [
    IdentifyTest(FileFormat.Fasta, 'simple.fasta'),
    IdentifyTest(FileFormat.Tabfile, 'simple.tsv'),
    IdentifyTest(FileFormat.Unknown, 'empty.txt'),
])
def test_write_sequences(test: IdentifyTest) -> None:
    test.validate()
