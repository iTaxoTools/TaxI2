from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

import pytest

from itaxotools.taxi2.pairs import SequencePairs

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class AcceptanceTest(NamedTuple):
    query: str
    reference: str
    output: str

    def validate(self, generated: SequencePairs):
        fixture_list = list(self.fixture())
        generated_list = list(generated)
        assert len(fixture_list) == len(generated_list)
        for pair in fixture_list:
            assert pair in generated_list


acceptance_tests = [
    AcceptanceTest("simple.query.tsv", "simple.reference.tsv", "simple.output.tsv")
]


@pytest.mark.xfail
@pytest.mark.parametrize("test", acceptance_tests)
def test_new_versus_reference(test: AcceptanceTest, tmp_path: Path) -> None:
    # query_path = TEST_DATA_DIR / test.query
    # reference_path = TEST_DATA_DIR / test.reference
    # fixed_path = TEST_DATA_DIR / test.output
    # output_path = tmp_path / test.output
    raise NotImplementedError()
