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
    input: tuple[str, str]
    solutions: list[tuple[str, str]]
    scores: tuple = (1, -1, -8, -1, -1, -1)

    def scores_from_tuple(self, scores: tuple[float]):
        return Scores(**{k: v for k, v in zip(Scores.defaults, scores)})

    def check(self, aligner_type: PairwiseAligner):
        scores = self.scores_from_tuple(self.scores)
        aligner = aligner_type(scores)
        x = Sequence('idx', self.input[0])
        y = Sequence('idy', self.input[1])
        ax, ay = aligner.align(SequencePair(x, y))
        print(ax, ay)
        assert ax.id == x.id
        assert ay.id == y.id
        assert len(ax.seq) == len(ay.seq)
        assert any(solution == (ax.seq, ay.seq) for solution in self.solutions)


aligner_types = [
    PairwiseAligner.Rust,
    PairwiseAligner.Biopython,
]


align_tests = [
    AlignTest(('TACTG', 'ACG'), [('TACTG', '-AC-G')], (1, 0, 0, 0, 0, 0)),
    AlignTest(('TACTG', 'ACG'), [('TACTG', '-ACG-')], (1, -1, -8, -1, -1, -1)),
    AlignTest(('TACTG', 'ACG'), [('TACTG', '-AC-G')], (1, -1, -1, -1, -1, -1)),
    AlignTest(('TACTG', 'ACG'), [('TACTG', '-ACG-')], (1, 0, -2, 0, 0, 0)),
    AlignTest(('TACTG', 'ACG'), [('TACTG', 'A-C-G')], (1, 0, 0, 0, -2, 0)),
    AlignTest(('TACTG', 'ACG'), [('TACTG', 'ACG--')], (0, 1, -1, 0, 0, 0)),

    AlignTest(('ATCG', 'ATAG'), [('ATC-G', 'AT-AG'), ('AT-CG', 'ATA-G'), ('ATCG', 'ATAG')], (1, 0, 0, 0, 0, 0)),
    AlignTest(('ATCG', 'ATAG'), [('ATC-G', 'AT-AG'), ('AT-CG', 'ATA-G')], (1, -1, 0, 0, 0, 0)),
    AlignTest(('ATCG', 'ATAG'), [('ATCG', 'ATAG')], (1, 0, -1, 0, 0, 0)),

    AlignTest(('ATCG', 'AG'), [('ATCG', 'A--G')], (1, 0, 0, 0, 0, 0)),
    AlignTest(('ATCG', 'AG'), [('ATCG', 'AG--'), ('ATCG', '--AG')], (1, 0, -2, 0, 0, 0)),
    AlignTest(('ATCG', 'AG'), [('ATCG', 'A--G')], (1, 0, -2, 0, -2, 0)),
    AlignTest(('ATCG', 'AG'), [('ATCG', 'A--G'), ('ATCG', '-AG-')], (1, 0, -2, 0, 0, -2)),
    AlignTest(('ATCG', 'AG'), [('ATCG', '-AG-')], (0, 0, -1, 0, 0, -1)),

    AlignTest(('ATATA', 'AAA'), [('ATATA', 'A-A-A')], (1, 0, 0, 0, 0, 0)),
    AlignTest(('ATATA', 'AAA'), [('ATATA', 'AAA--'), ('ATATA', '--AAA')], (1, 0, -1, 0, 0, 0)),

    AlignTest(('ATATATATATA', 'ATTA'), [('ATATATATATA', 'AT-------TA')], (1, 0, 0, 0.1, 0, 0)),
]


@pytest.mark.parametrize("aligner_type", aligner_types)
@pytest.mark.parametrize("test", align_tests)
def test_align(aligner_type: PairwiseAligner, test: AlignTest) -> None:
    test.check(aligner_type)
