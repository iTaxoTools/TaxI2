from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import NamedTuple

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .types import Container, Type
from .handlers import FileHandler


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
    def _iter_read(self) -> iter[Sequence]:
        with open(self.path, 'r') as handle:
            yield  # ready
            for data in SimpleFastaParser(handle):
                yield Sequence(*data)

    def _iter_write(self) -> iter[Sequence]:
        raise NotImplementedError()


class Genbank(SequenceHandler):
    def _iter_read(self) -> iter[Sequence]:
        # Bio.GenBank.Scanner
        file = SeqIO.parse(self.path, 'genbank')
        yield  # ready
        for data in file:
            yield Sequence(data.id, data.seq)

    def _iter_write(self) -> iter[Sequence]:
        raise NotImplementedError()


class Tabular(SequenceHandler):
    subhandler = FileHandler.Tabular

    def _open_readable(
        self,
        idHeader: str = None,
        seqHeader: str = None,
        hasHeader: bool = False,
        idColumn: int = 0,
        seqColumn: int = 1,
    ):
        if idHeader and seqHeader:
            columns = (idHeader, seqHeader)
            hasHeader = True
        else:
            columns = (idColumn, seqColumn)

        self.columns = columns
        self.hasHeader = hasHeader
        super()._open_readable()

    def _iter_read(self) -> iter[Sequence]:
        with self.subhandler(
            self.path,
            has_headers=self.hasHeader,
            columns=self.columns,
            get_all_columns=True,
        ) as rows:
            headers = rows.headers
            extras = dict()
            yield
            for row in rows:
                id = row[0]
                seq = row[1]
                if headers is not None:
                    extras = { k: v for (k, v) in zip(headers[2:], row[2:]) }
                yield Sequence(id, seq, extras)


class Tabfile(SequenceHandler.Tabular, SequenceHandler):
    subhandler = FileHandler.Tabular.Tabfile

    def _open_writable(
        self,
        idHeader: str = None,
        seqHeader: str = None,
        hasHeader: bool = False,
    ):
        if idHeader and seqHeader:
            hasHeader = True
        self.idHeader = idHeader
        self.seqHeader = seqHeader
        self.hasHeader = hasHeader
        super()._open_writable()

    def _iter_write(self, idHeader='seqid', seqHeader='sequence'):
        with self.subhandler(self.path, 'w') as file:
            try:
                sequence = yield
                if self.hasHeader:
                    extraHeaders = tuple(sequence.extras.keys())
                    file.write((self.idHeader,) + extraHeaders + (self.seqHeader,))
                while True:
                    extras = tuple(sequence.extras.values())
                    file.write((sequence.id,) + extras + (sequence.seq,))
                    sequence = yield
            except GeneratorExit:
                return


class Excel(SequenceHandler.Tabular, SequenceHandler):
    subhandler = FileHandler.Tabular.Excel

    def _iter_write(self) -> iter[Sequence]:
        raise NotImplementedError()
