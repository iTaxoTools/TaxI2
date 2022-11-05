from __future__ import annotations

from pathlib import Path
from typing import NamedTuple
from contextlib import contextmanager
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from openpyxl import load_workbook

from .types import Container, Type


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

    def __init__(self, path: Path):
        self.path = path

    def read(self, *args, **kwargs) -> iter[Sequence]:
        raise NotImplementedError()

    def write(self, distances: iter[Sequence], *args, **kwargs) -> None:
        raise NotImplementedError()


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
    def iter_rows(self):
        raise NotImplementedError()

    @contextmanager
    def open(self):
        yield self.iter_rows()

    def read(
        self,
        idHeader: str = None,
        seqHeader: str = None,
        hasHeader: bool = False,
        idColumn: int = 0,
        seqColumn: int = 1,
    ) -> iter[Sequence]:

        with self.open() as rows:

            # Checking headers
            if idHeader and seqHeader:
                hasHeader = True
            if hasHeader:
                headers = next(rows)
                idColumn, seqColumn = headers.index(idHeader), headers.index(seqHeader)

            # Getting id and seq
            for row in rows:
                if len(row) <= 1:
                    continue
                id = row[idColumn]
                seq = row[seqColumn]

                extras = {}
                if hasHeader:
                    extras = {
                        k: v for c, (k, v) in enumerate(zip(headers, row))
                        if c not in (idColumn, seqColumn)
                    }

                yield Sequence(id, seq, extras)


class Tabfile(Tabular, SequenceFile):
    def iter_rows(self) -> iter[tuple[str, ...]]:
        with open(self.path) as file:
            for line in file:
                yield line.strip().split('\t')


class Excel(Tabular, SequenceFile):
    def iter_rows(self) -> iter[tuple[str, ...]]:
        wb = load_workbook(filename=self.path, read_only=True)
        ws = wb.worksheets[0]
        for row in ws.iter_rows(values_only=True):
            row = list(row)
            while row and row[-1] is None:
                del row[-1]
            yield [x if x else '' for x in row]
        wb.close
