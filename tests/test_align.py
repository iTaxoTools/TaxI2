from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest
from utility import assert_eq_files

from itaxotools.taxi3.align import PairwiseAligner, Scores
from itaxotools.taxi3.pairs import (
    SequencePair, SequencePairFile, SequencePairs)
from itaxotools.taxi3.sequences import Sequence, Sequences

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class AlignTest(NamedTuple):
    seq_x: str
    seq_y: str
    aligned_x: str
    aligned_y: str
    scores: dict = {}

    def check(self, aligner_type: PairwiseAligner):
        scores = Scores(**self.scores)
        aligner = aligner_type(scores)
        x = Sequence('idx', self.seq_x)
        y = Sequence('idy', self.seq_y)
        ax, ay = aligner.align(SequencePair(x, y))
        print(ax, ay)
        assert ax.id == x.id
        assert ay.id == y.id
        assert len(ax.seq) == len(ay.seq)
        assert ax.seq == self.aligned_x
        assert ay.seq == self.aligned_y


aligner_types = [
    PairwiseAligner.Rust,
    PairwiseAligner.Biopython,
]


align_tests = [
    AlignTest('TACTG', 'ACG', 'TACTG', '-ACG-', {}),
    AlignTest('TACTG', 'ACG', 'TACTG', '-AC-G', dict(gap_penalty=-1)),
    AlignTest('TACTG', 'ACG----', 'TACTG---', '-ACG----', {}),
    AlignTest('---TACTG', 'ACG-----', '------TACTG', 'ACG--------', {}),
    AlignTest('TACTG---', '-----ACG', 'TACTG------', '--------ACG', {}),
    AlignTest('ACTG--', '--TGAC', 'ACTG--', '--TGAC', {}),
]


@pytest.mark.parametrize("aligner_type", aligner_types)
@pytest.mark.parametrize("test", align_tests)
def test_align(aligner_type: PairwiseAligner, test: AlignTest) -> None:
    test.check(aligner_type)
