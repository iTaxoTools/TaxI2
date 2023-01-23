from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import NamedTuple

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .types import Container, Type
from .handlers import FileHandler, ReadHandle, WriteHandle
from .encoding import sanitize


class Sequence(NamedTuple):
    id: str
    seq: str
    extras: dict[str, str] = dict()

    _tr_normalize = str.maketrans('?', 'N', '-')

    def normalize(self):
        return Sequence(self.id, self.seq.translate(self._tr_normalize).upper(), self.extras)


class Sequences(Container[Sequence]):
    @classmethod
    def fromPath(cls, path: Path, handler: SequenceHandler, *args, **kwargs) -> Sequences:
        return cls(handler, path, *args, **kwargs)

    def normalize(self) -> Sequences:
        return Sequences(lambda: (seq.normalize() for seq in self))


class SequenceHandler(FileHandler[Sequence]):
    pass


class Fasta(SequenceHandler):
    def _iter_read(self) -> ReadHandle[Sequence]:
        with open(self.path, 'r') as handle:
            yield self
            for data in SimpleFastaParser(handle):
                yield Sequence(*data)

    def _iter_write(self) -> WriteHandle[Sequence]:
        raise NotImplementedError()


class Genbank(SequenceHandler):
    def _iter_read(self) -> ReadHandle[Sequence]:
        # Bio.GenBank.Scanner
        file = SeqIO.parse(self.path, 'genbank')
        yield self
        for data in file:
            yield Sequence(data.id, data.seq)

    def _iter_write(self) -> WriteHandle[Sequence]:
        raise NotImplementedError()


class Tabular(SequenceHandler):
    subhandler = FileHandler.Tabular

    def _iter_read(
        self,
        idHeader: str = None,
        seqHeader: str = None,
        hasHeader: bool = False,
        idColumn: int = 0,
        seqColumn: int = 1,
    ) -> ReadHandle[Sequence]:

        if idHeader and seqHeader:
            columns = (idHeader, seqHeader)
            hasHeader = True
        else:
            columns = (idColumn, seqColumn)

        with self.subhandler(
            self.path,
            has_headers=hasHeader,
            columns=columns,
            get_all_columns=True,
        ) as rows:

            headers = rows.headers
            if headers is not None:
                headers = [sanitize(header) for header in headers]
            extras = dict()
            yield self
            for row in rows:
                id = sanitize(row[0])
                seq = row[1]
                if headers is not None:
                    extras = { k: v for (k, v) in zip(headers[2:], row[2:]) }
                yield Sequence(id, seq, extras)


class Tabfile(SequenceHandler.Tabular, SequenceHandler):
    subhandler = FileHandler.Tabular.Tabfile

    def _iter_write(
        self,
        idHeader: str = None,
        seqHeader: str = None,
        hasHeader: bool = False,
    ) -> WriteHandle[Sequence]:

        wrote_headers = False
        if idHeader and seqHeader:
            hasHeader = True

        with self.subhandler(self.path, 'w') as file:
            try:
                sequence = yield
                if hasHeader:
                    extraHeaders = tuple(sequence.extras.keys())
                    file.write((idHeader,) + extraHeaders + (seqHeader,))
                    wrote_headers = True
                while True:
                    extras = tuple(sequence.extras.values())
                    file.write((sequence.id,) + extras + (sequence.seq,))
                    sequence = yield
            except GeneratorExit:
                if hasHeader and not wrote_headers:
                    file.write((idHeader, seqHeader))


class Excel(SequenceHandler.Tabular, SequenceHandler):
    subhandler = FileHandler.Tabular.Excel

    def _iter_write(self) -> WriteHandle[Sequence]:
        raise NotImplementedError()
