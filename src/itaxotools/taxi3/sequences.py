from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import NamedTuple

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .types import Container, Type
from .handlers import FileHandler, TabularHandler, TabfileHandler, ExcelHandler


class Sequence(NamedTuple):
    id: str
    seq: str
    extras: dict[str, str] = dict()

    _tr_normalize = str.maketrans('?', 'N', '-')

    def normalize(self):
        return Sequence(self.id, self.seq.translate(self._tr_normalize).upper(), self.extras)


class Sequences(Container[Sequence]):
    @classmethod
    def fromFile(cls, file: SequenceFile, *args, **kwargs) -> Sequences:
        return cls(file.open, *args, **kwargs)

    def normalize(self) -> Sequences:
        return Sequences(lambda: (seq.normalize() for seq in self))


class Handler(Type, FileHandler[Sequence]):
    pass


class Fasta(Handler):
    def _iter_read(self) -> iter[Sequence]:
        with open(self.path, 'r') as handle:
            yield  # ready
            for data in SimpleFastaParser(handle):
                yield Sequence(*data)


class Genbank(Handler):
    def _iter_read(self) -> iter[Sequence]:
        # Bio.GenBank.Scanner
        file = SeqIO.parse(self.path, 'genbank')
        yield  # ready
        for data in file:
            yield Sequence(data.id, data.seq)


class Tabular(Handler):
    subhandler = TabularHandler

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


class Tabfile(Handler.Tabular, Handler):
    subhandler = TabfileHandler

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


class Excel(Handler.Tabular, Handler):
    subhandler = ExcelHandler


class SequenceFile(Type):
    handler = Handler

    def __init__(self, path: Path):
        self.path = path

    def open(self, *args, **kwargs) -> Handler:
        return self.handler(self.path, *args, **kwargs)

    @classmethod
    def identifyFile(cls, file):
        file_map = {'tab': Tabfile, 'csv': Excel, 'fas': Fasta, 'gbk': Genbank}
        file_type = str(file).strip().split('.')[1]
        return file_map[file_type](file)


class Fasta(SequenceFile):
    handler = Handler.Fasta


class Genbank(SequenceFile):
    handler = Handler.Genbank


class Tabfile(SequenceFile):
    handler = Handler.Tabfile


class Excel(SequenceFile):
    handler = Handler.Excel
