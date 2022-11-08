from __future__ import annotations

import multiprocessing
from warnings import warn

from Bio.Align import PairwiseAligner as BioPairwiseAligner
from Bio.Seq import reverse_complement

from itaxotools.taxi3.library import calculate_distances as calc

from .pairs import SequencePair, SequencePairs
from .sequences import Sequence
from .types import Type


class Scores(dict[str, int]):
    """Can access keys like attributes"""

    defaults = dict(
        match_score = 1,
        mismatch_score = -1,
        internal_open_gap_score = -8,
        internal_extend_gap_score = -1,
        end_open_gap_score = -1,
        end_extend_gap_score = -1,
    )

    def __init__(self, **kwargs):
        super().__init__(self.defaults | kwargs)
        self.__dict__ = self

    def __repr__(self):
        attrs = ', '.join(f'{k}={v}' for k, v in self.items())
        return f'<{type(self).__name__}: {attrs}>'


class PairwiseAligner(Type):
    def __init__(self, scores: Scores = None):
        self.scores = scores or Scores()

    def align(self, pair: SequencePair) -> SequencePair:
        raise NotImplementedError()

    def align_pairs_parallel(self, pairs: SequencePairs) -> SequencePairs:
        with multiprocessing.Pool(processes=4, maxtasksperchild=10) as pool:
            for x in pool.imap(self.align, pairs, chunksize=1000):
                yield x

    def align_pairs(self, pairs: SequencePairs) -> SequencePairs:
        return SequencePairs((self.align(pair) for pair in pairs))


class Rust(PairwiseAligner):

    def __init__(self, scores: Scores = None):
        super().__init__(scores)
        warn('PairwiseAligner.Rust does not always find the best alignment!', RuntimeWarning)
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

    def _format_pretty(self, alignment):
        # Adjusted from Bio.Align.PairwiseAlignment._format_pretty
        seq1 = alignment._convert_sequence_string(alignment.target)
        if seq1 is None:
            return alignment._format_generalized()
        seq2 = alignment._convert_sequence_string(alignment.query)
        if seq2 is None:
            return alignment._format_generalized()
        n2 = len(seq2)
        aligned_seq1 = ""
        aligned_seq2 = ""
        pattern = ""
        path = alignment.path
        if path[0][1] > path[-1][1]:  # mapped to reverse strand
            path = tuple((c1, n2 - c2) for (c1, c2) in path)
            seq2 = reverse_complement(seq2)
        end1, end2 = path[0]
        if end1 > 0 or end2 > 0:
            end = max(end1, end2)
            aligned_seq1 += " " * (end - end1) + seq1[:end1]
            aligned_seq2 += " " * (end - end2) + seq2[:end2]
            pattern += " " * end
        start1 = end1
        start2 = end2
        for end1, end2 in path[1:]:
            if end1 == start1:
                gap = end2 - start2
                aligned_seq1 += "-" * gap
                aligned_seq2 += seq2[start2:end2]
                pattern += "-" * gap
            elif end2 == start2:
                gap = end1 - start1
                aligned_seq1 += seq1[start1:end1]
                aligned_seq2 += "-" * gap
                pattern += "-" * gap
            else:
                s1 = seq1[start1:end1]
                s2 = seq2[start2:end2]
                aligned_seq1 += s1
                aligned_seq2 += s2
                for c1, c2 in zip(s1, s2):
                    if c1 == c2:
                        pattern += "|"
                    else:
                        pattern += "."
            start1 = end1
            start2 = end2
        aligned_seq1 += seq1[end1:]
        aligned_seq2 += seq2[end2:]
        return (aligned_seq1, pattern, aligned_seq2)

    def align(self, pair: SequencePair) -> SequencePair:
        alignments = self.aligner.align(pair.x.seq, pair.y.seq)
        aligned_x, _, aligned_y = self._format_pretty(alignments[0])
        return SequencePair(
            Sequence(pair.x.id, aligned_x, pair.x.extras),
            Sequence(pair.y.id, aligned_y, pair.y.extras),
        )
