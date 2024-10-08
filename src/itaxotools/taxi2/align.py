from __future__ import annotations

import multiprocessing
from typing import Iterator
from warnings import warn

from Bio.Align import Alignment
from Bio.Align import PairwiseAligner as BioPairwiseAligner

from itaxotools import calculate_distances as calc

from .pairs import SequencePair, SequencePairs
from .sequences import Sequence
from .types import Type


class Scores(dict[str, int]):
    """Can access keys like attributes"""

    defaults = dict(
        match_score=1,
        mismatch_score=-1,
        internal_open_gap_score=-8,
        internal_extend_gap_score=-1,
        end_open_gap_score=-1,
        end_extend_gap_score=-1,
    )

    def __init__(self, **kwargs):
        super().__init__(self.defaults | kwargs)
        self.__dict__ = self

    def __repr__(self):
        attrs = ", ".join(f"{k}={v}" for k, v in self.items())
        return f"<{type(self).__name__}: {attrs}>"


class PairwiseAligner(Type):
    def __init__(self, scores: Scores = None):
        self.scores = scores or Scores()

    def align(self, pair: SequencePair) -> SequencePair:
        raise NotImplementedError()

    def align_pairs_parallel(self, pairs: SequencePairs) -> Iterator[SequencePair]:
        with multiprocessing.Pool(processes=4, maxtasksperchild=10) as pool:
            for x in pool.imap(self.align, pairs, chunksize=1000):
                yield x

    def align_pairs(self, pairs: SequencePairs) -> SequencePairs:
        return SequencePairs((self.align(pair) for pair in pairs))


class Rust(PairwiseAligner):
    def __init__(self, scores: Scores = None):
        super().__init__(scores)
        warn(
            "PairwiseAligner.Rust does not always find the best alignment!",
            RuntimeWarning,
        )
        self.aligner = calc.make_aligner(**self.scores)

    def align(self, pair: SequencePair) -> SequencePair:
        alignments = calc.align_seq(self.aligner, pair.x.seq, pair.y.seq)
        aligned_x, aligned_y = alignments
        return SequencePair(
            Sequence(pair.x.id, aligned_x, pair.x.extras),
            Sequence(pair.y.id, aligned_y, pair.y.extras),
        )


class Biopython(PairwiseAligner):
    def __init__(self, scores: Scores = None):
        super().__init__(scores)
        self.aligner = BioPairwiseAligner(**self.scores)

    def _format_pretty(self, alignment: Alignment):
        # Adjusted from Bio.Align.Alignment._format_generalized
        seq1, seq2 = alignment.sequences
        aligned_seq1 = []
        aligned_seq2 = []
        pattern = []
        end1, end2 = alignment.coordinates[:, 0]
        if end1 > 0 or end2 > 0:
            if end1 <= end2:
                for c2 in seq2[: end2 - end1]:
                    s2 = str(c2)
                    s1 = " " * len(s2)
                    aligned_seq1.append(s1)
                    aligned_seq2.append(s2)
                    pattern.append(s1)
            else:  # end1 > end2
                for c1 in seq1[: end1 - end2]:
                    s1 = str(c1)
                    s2 = " " * len(s1)
                    aligned_seq1.append(s1)
                    aligned_seq2.append(s2)
                    pattern.append(s2)
        start1 = end1
        start2 = end2
        for end1, end2 in alignment.coordinates[:, 1:].transpose():
            if end1 == start1:
                for c2 in seq2[start2:end2]:
                    s2 = str(c2)
                    s1 = "-" * len(s2)
                    aligned_seq1.append(s1)
                    aligned_seq2.append(s2)
                    pattern.append(s1)
                start2 = end2
            elif end2 == start2:
                for c1 in seq1[start1:end1]:
                    s1 = str(c1)
                    s2 = "-" * len(s1)
                    aligned_seq1.append(s1)
                    aligned_seq2.append(s2)
                    pattern.append(s2)
                start1 = end1
            else:
                t1 = seq1[start1:end1]
                t2 = seq2[start2:end2]
                if len(t1) != len(t2):
                    raise ValueError("Unequal step sizes in alignment")
                for c1, c2 in zip(t1, t2):
                    s1 = str(c1)
                    s2 = str(c2)
                    m1 = len(s1)
                    m2 = len(s2)
                    if c1 == c2:
                        p = "|"
                    else:
                        p = "."
                    if m1 < m2:
                        space = (m2 - m1) * " "
                        s1 += space
                        pattern.append(p * m1 + space)
                    elif m1 > m2:
                        space = (m1 - m2) * " "
                        s2 += space
                        pattern.append(p * m2 + space)
                    else:
                        pattern.append(p * m1)
                    aligned_seq1.append(s1)
                    aligned_seq2.append(s2)
                start1 = end1
                start2 = end2
        aligned_seq1 = "".join(aligned_seq1)
        aligned_seq2 = "".join(aligned_seq2)
        pattern = "".join(pattern)
        return (aligned_seq1, pattern, aligned_seq2)

    def align(self, pair: SequencePair) -> SequencePair:
        alignments = self.aligner.align(pair.x.seq, pair.y.seq)
        aligned_x, _, aligned_y = self._format_pretty(alignments[0])
        return SequencePair(
            Sequence(pair.x.id, aligned_x, pair.x.extras),
            Sequence(pair.y.id, aligned_y, pair.y.extras),
        )
