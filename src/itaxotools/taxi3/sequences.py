from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import NamedTuple

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .types import Container, Type
from .tabular import Tabular as TabularProtocol


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
        return cls(file.read, *args, **kwargs)

    def normalize(self) -> Sequences:
        return Sequences(lambda: (seq.normalize() for seq in self))


class SequenceFile(Type):
    """Handlers for sequence files"""

    def __init__(self, path: Path) -> object:
        self.path = path

    def read(self, *args, **kwargs) -> iter[Sequence]:
        raise NotImplementedError()

    def write(self, distances: iter[Sequence], *args, **kwargs) -> None:
        raise NotImplementedError()

    @classmethod
    def identifyFile(cls, file):
        file_map = {'tab': Tabfile, 'csv': Excel, 'fas': Fasta, 'gbk': Genbank}
        file_type = str(file).strip().split('.')[1]
        return file_map[file_type](file)


class SequenceFileHandler:
    def __init__(self, file: SequenceFile):
        self.coroutine = file.iter_write()
        next(self.coroutine)

    def write(self, sequence: Sequence) -> None:
        self.coroutine.send(sequence)

    def close(self) -> None:
        self.coroutine.close()


class Fasta(SequenceFile):
    def read(self) -> iter[Sequence]:
        with open(self.path, 'r') as handle:
            for data in SimpleFastaParser(handle):
                yield Sequence(*data)


class Genbank(SequenceFile):
    def read(self) -> iter[Sequence]:
        # Bio.GenBank.Scanner
        for data in SeqIO.parse(self.path, 'genbank'):
            yield Sequence(data.id, data.seq)


class Tabular(SequenceFile):
    protocol = TabularProtocol

    @contextmanager
    def open(self, mode='r'):
        if mode == 'r':
            yield self.iter_rows()
        elif mode == 'w':
            handler = SequenceFileHandler(self)
            try:
                yield handler
            finally:
                handler.close()
        else:
            raise Exception('Unknown mode')

    def read(
        self,
        idHeader: str = None,
        seqHeader: str = None,
        hasHeader: bool = False,
        idColumn: int = 0,
        seqColumn: int = 1,
    ) -> iter[Sequence]:

        if idHeader and seqHeader:
            columns = (idHeader, seqHeader)
            hasHeader = True
        else:
            columns = (idColumn, seqColumn)

        with self.protocol.open(self.path, has_headers=hasHeader, columns=columns, get_all_columns=True) as rows:
            headers = rows.headers
            extras = dict()
            for row in rows:
                id = row[0]
                seq = row[1]
                if headers is not None:
                    extras = { k: v for (k, v) in zip(headers[2:], row[2:]) }
                yield Sequence(id, seq, extras)

    def getHeader(self):
        return cls.protocol.headers(self.path)


class Tabfile(Tabular, SequenceFile):
    protocol = TabularProtocol.Tabfile

    def iter_write(self, idHeader='seqid', seqHeader='sequence'):
        with open(self.path, 'w') as file:
            try:
                sequence = yield
                extraHeaders = sequence.extras.keys()
                extraHeadersString = '\t'.join(extraHeaders)
                file.write(f'{idHeader}\t{extraHeadersString}\t{seqHeader}\n')
                while True:
                    extras = sequence.extras.values()
                    extrasString = '\t'.join(extras)
                    file.write(f'{sequence.id}\t{extrasString}\t{sequence.seq}\n')
                    sequence = yield  # Yield expression
            except GeneratorExit:
                return


class Excel(Tabular, SequenceFile):
    protocol = TabularProtocol.Excel
