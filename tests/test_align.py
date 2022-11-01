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
        print(ax, ay, scores)
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
    AlignTest(('ATCG', 'AG'), [('ATCG', '-AG-')], (0, 0, -1, 0, 0, -1)),

    AlignTest(('ATATA', 'AAA'), [('ATATA', 'A-A-A')], (1, 0, 0, 0, 0, 0)),
    AlignTest(('ATATA', 'AAA'), [('ATATA', 'AAA--'), ('ATATA', '--AAA')], (1, 0, -1, 0, 0, 0)),
    #AlignTest(('ATATATATATA', 'ATTA'), [('ATATATATATA', 'AT-------TA')], (10, 0, 0, 1, 0, 0)),

    # simple match
    AlignTest(('ATCG', 'ATCG'), [('ATCG', 'ATCG')], (1, 0, 0, 0, 0, 0)),
    AlignTest(('ATCG', 'AT'), [('ATCG', 'AT--')], (1, 0, 0, 0, 0, 0)),
    AlignTest(('ATCG', 'CG'), [('ATCG', '--CG')], (1, 0, 0, 0, 0, 0)),
    AlignTest(('ATCG', 'TC'), [('ATCG', '-TC-')], (1, 0, 0, 0, 0, 0)),

    # mismatch score
    AlignTest(('ATCG', 'GCTA'), [('ATCG', 'GCTA')], (1, 1, 0, 0, 0, 0)),
    AlignTest(('ATCG', 'ATCG'), [('ATCG-', '-ATCG'), ('-ATCG', 'ATCG-')], (0, 1, 0, 0, 0, 0)),
    AlignTest(('ATC', 'AGC'), [('AT-C', 'A-GC'), ('A-TC', 'AG-C')], (1, -1, 0, 0, 0, 0)),
    AlignTest(('AAT', 'AAC'), [('AAT-', 'AA-C'), ('AA-T', 'AAC-')], (1, -1, 0, 0, 0, 0)),
    AlignTest(('TAA', 'CAA'), [('-TAA', 'C-AA'), ('T-AA', '-CAA')], (1, -1, 0, 0, 0, 0)),

    # gap penalty: open internal
    AlignTest(('AAT', 'AAC'), [('AAT', 'AAC')], (1, 0, -1, 0, 0, 0)),
    AlignTest(('TAA', 'CAA'), [('TAA', 'CAA')], (1, 0, -1, 0, 0, 0)),
    AlignTest(('ATC', 'AGC'), [('ATC', 'AGC')], (1, 0, -1, 0, 0, 0)),
    AlignTest(('ATC', 'AGC'), [('ATC', 'AGC')], (1, -1, -1, 0, 0, 0)),
    #AlignTest(('AAATTTAAA', 'AAACCCAAA'), [('AAA---TTTAAA', 'AAACCC---AAA'), ('AAATTT---AAA', 'AAA---CCCAAA')], (1, -1, -1, 0, 0, 0)),
    AlignTest(('AAATTTAAA', 'AAACCCAAA'), [('AAA---TTTAAA', 'AAACCC---AAA'), ('AAATTT---AAA', 'AAA---CCCAAA')], (1, -2, -1, 0, 0, 0)),
    AlignTest(('AAATTTAAA', 'AAACCCAAA'), [('AAATTTAAA', 'AAACCCAAA'), ('------AAATTTAAA', 'AAACCCAAA------'), ('AAATTTAAA------', '------AAACCCAAA')], (1, -1, -2, 0, 0, 0)),
    AlignTest(('AAACTAAA', 'AAATGAAA'), [('AAACT-AAA', 'AAA-TGAAA')], (1, -1, -1, 0, 0, 0)),
    AlignTest(('AAACTAAA', 'AAATGAAA'), [('AAACTAAA', 'AAATGAAA')], (1, -1, -2, 0, 0, 0)),

    # gap penalty: end & internal
    AlignTest(('AAA', 'TTT'), [('AAA', 'TTT')], (1, 0, -1, 0, -1, 0)),

    # gap penalty: extend internal
    AlignTest(('ATACCGG', 'ATAGG'), [('ATACCGG', 'ATA--GG')], (1, -1, 0, 0, 0, 0)),
    AlignTest(('ATACCGG', 'ATAGG'), [('ATAC-CGG', 'ATA-G-G-')], (1, -1, 0, -2, 0, 0)),

]


@pytest.mark.parametrize("aligner_type", aligner_types)
@pytest.mark.parametrize("test", align_tests)
def test_align(aligner_type: PairwiseAligner, test: AlignTest) -> None:
    test.check(aligner_type)
