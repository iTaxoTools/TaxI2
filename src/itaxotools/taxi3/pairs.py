from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

from .sequences import Sequence, Sequences
from .types import Container, Type


class SequencePair(NamedTuple):
    x: Sequence
    y: Sequence


class SequencePairs(Container[SequencePair]):
    @classmethod
    def fromFile(cls, file: SequencePairFile, *args, **kwargs) -> SequencePairs:
        return cls(file.read, *args, **kwargs)

    @classmethod
    def fromProduct(cls, xs: Sequences, ys: Sequences) -> SequencePairs:
        return cls(SequencePair(x, y) for x in xs for y in ys)


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
        with open(self.path, 'r') as f:
            for line in f:
                lineData = line[:-1].split('\t')
                print(lineData)
                idx, idy, seqX, seqY = lineData[0], lineData[1], lineData[2], lineData[3]
                yield SequencePair(Sequence(idx, seqX),Sequence(idy, seqY))

    def write(self, pairs: iter[SequencePair], *args, **kwargs) -> None:
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
        with open(self.path, 'r') as f:
            buffer = []
            for line in f:
                line = line.strip()
                # Skipping blank lines
                if not line:
                    continue

                buffer.append(line)
                if len(buffer) == 4:
                    idx, idy = buffer[0].strip().split('/')
                    seqX, seqY = buffer[1], buffer[3] # Skipping pos 2 as it is the format
                    buffer = []
                    yield SequencePair(Sequence(idx.strip(), seqX.strip()),Sequence(idy.strip(), seqY.strip()))


    def write(self, pairs: iter[SequencePair], *args, **kwargs) -> None:
        raise NotImplementedError()
