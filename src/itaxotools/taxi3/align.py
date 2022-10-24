from __future__ import annotations

from typing import Callable, Iterable, NamedTuple
from enum import Enum, auto

from .sequences import Sequence
from .types import Type, Container


class SequencePair(NamedTuple):
    x: Sequence
    y: Sequence


class SequencePairs(Container[SequencePair]):
    @classmethod
    def fromFile(cls, file: SequencePairFile, *args, **kwargs) -> SequencePairs:
        return cls(file.read, *args, **kwargs)


class Scores(dict):
    """Can access keys like attributes"""

    defaults = dict(
        gap_penalty = -8,
        gap_extend_penalty = -1,
        end_gap_penalty = -1,
        end_gap_extend_penalty = -1,
        match_score = 1,
        mismatch_score = -1,
    )

    def __init__(self, **kwargs):
        super().__init__(self.defaults | kwargs)
        self.__dict__ = self


class PairwiseAligner:
    def __init__(self, scores: Scores = None):
        self.scores = scores or Scores()

    def align(self, x: Sequence, y: Sequence) -> SequencePair:
        raise NotImplementedError()


class SequencePairFile(Type):
    """Handlers for pairwise-aligned sequence files"""

    def __init__(self, path: Path):
        self.path = path

    def read(self, *args, **kwargs) -> iter[SequencePair]:
        raise NotImplementedError()

    def write(self, pairs: iter[SequencePair], *args, **kwargs) -> None:
        raise NotImplementedError()


class Tabfile(SequencePairFile):

    def read(self) -> iter[SequencePair]:
        raise NotImplementedError()

    def write(self, pairs: iter[SequencePair]) -> None:
        raise NotImplementedError()


class Formatted(SequencePairFile):

    @staticmethod
    def _format_char(x: str, y: str) -> str:
        if x == y:
            return '|'
        if x == '-' or y == '-':
            return ' '
        return '.'

    @classmethod
    def _format(self, x: str, y: str) -> str:
        return ''.join((self._format_char(a, b) for a, b in zip(x, y)))

    def read(self) -> iter[SequencePair]:
        raise NotImplementedError()

    def write(self, pairs: iter[SequencePair]) -> None:
        raise NotImplementedError()
