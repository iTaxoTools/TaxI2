from __future__ import annotations

from pathlib import Path
from typing import NamedTuple, TextIO

from .handlers import FileHandler, ReadHandle, WriteHandle
from .sequences import Sequence, Sequences
from .types import Container


class SequencePair(NamedTuple):
    x: Sequence
    y: Sequence


class SequencePairs(Container[SequencePair]):
    @classmethod
    def fromPath(
        cls, path: Path, handler: SequencePairHandler, *args, **kwargs
    ) -> SequencePairs:
        return cls(handler, path, *args, **kwargs)

    @classmethod
    def fromProduct(cls, xs: Sequences, ys: Sequences) -> SequencePairs:
        return cls(lambda: (SequencePair(x, y) for x in xs for y in ys))


class SequencePairHandler(FileHandler[SequencePair]):
    pass


class Tabfile(SequencePairHandler):
    def _iter_read(self) -> ReadHandle[SequencePair]:
        with FileHandler.Tabfile(self.path, "r", has_headers=True) as file:
            yield self
            for idx, idy, seqX, seqY in file:
                yield SequencePair(Sequence(idx, seqX), Sequence(idy, seqY))

    def _iter_write(self) -> WriteHandle[SequencePair]:
        with FileHandler.Tabfile(
            self.path, "w", columns=["idx", "idy", "seqx", "seqy"]
        ) as file:
            try:
                while True:
                    pair = yield
                    file.write((pair.x.id, pair.y.id, pair.x.seq, pair.y.seq))
            except GeneratorExit:
                return


class Formatted(SequencePairHandler):
    @staticmethod
    def _format_char(x: str, y: str) -> str:
        if x == y and x != "-" and y != "-":
            return "|"
        if x == "-" or y == "-":
            return "-"
        return "."

    @classmethod
    def _format(self, x: str, y: str) -> str:
        return "".join((self._format_char(a, b) for a, b in zip(x, y)))

    def _iter_read(self) -> ReadHandle[SequencePair]:
        with open(self.path, "r") as file:
            yield self
            while lines := self._read_lines(file):
                idx, idy = lines[0].split(" / ")
                seqX, seqY = lines[1], lines[3]
                yield SequencePair(Sequence(idx, seqX), Sequence(idy, seqY))

    def _read_lines(self, file: TextIO):
        lines = []
        for _ in range(5):
            line = file.readline().strip()
            lines.append(line)
        if not any(lines):
            return []
        return lines

    def _iter_write(self) -> WriteHandle[SequencePair]:
        with open(self.path, "w") as file:
            try:
                pair = yield
                self._write_lines(file, pair)
                while True:
                    pair = yield
                    file.write("\n")
                    self._write_lines(file, pair)
            except GeneratorExit:
                return

    def _write_lines(self, file: TextIO, pair: SequencePair):
        file.write(f"{pair.x.id} / {pair.y.id}\n")
        file.write(f"{pair.x.seq}\n")
        file.write(f"{self._format(pair.x.seq, pair.y.seq)}\n")
        file.write(f"{pair.y.seq}\n")
