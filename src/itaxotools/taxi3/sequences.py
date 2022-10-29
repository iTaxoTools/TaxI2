from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from openpyxl import load_workbook

from .types import Container, Type


class Sequence(NamedTuple):
    id: str
    seq: str


class Sequences(Container[Sequence]):
    @classmethod
    def fromFile(cls, file: SequenceFile, *args, **kwargs) -> Sequences:
        return cls(file.read, *args, **kwargs)


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


class Tabfile(SequenceFile):
    def read(
        self,
        idHeader: str = None,
        seqHeader: str = None,
        hasHeader: bool = False,
        idColumn: int = 0,
        seqColumn: int = 1,
    ) -> iter[Sequence]:
        with open(self.path) as f:
            # Checking headers
            if idHeader and seqHeader:
                hasHeader = True
            if hasHeader:
                headerLine = f.readline().strip().split('\t', maxsplit=-1)
                idColumn, seqColumn = headerLine.index(idHeader), headerLine.index(seqHeader)

            # Getting id and seq
            for line in f:
                if len(line) <= 1:
                    continue
                data = line.strip().split('\t', maxsplit=-1)
                yield Sequence(data[idColumn], data[seqColumn])


class Excel(SequenceFile):
    def read(
        self,
        id: str = None,
        seq: str = None,
    ) -> iter[Sequence]:
        wb = load_workbook(filename=self.path, read_only=True)
        ws = wb['Sheet 1']
        idColumn = None
        seqColumn = None
        numCells = 0
        numRows = 0
        for row in ws.rows:
            lenRow = 0
            for cell in row:
                if cell.value is not None:
                    if cell.value == 'seqid':
                        idColumn = cell.column
                    elif cell.value == 'sequences':
                        seqColumn = cell.column
                    lenRow += 1
            if lenRow > 0:
                numRows = lenRow
            numCells += 1

        if numRows == 2:
            for x in range(1, numCells):
                id = ws.cell(row=x, column=1).value
                seq = ws.cell(row=x, column=2).value
                if id and seq:
                    yield Sequence(id, seq)

        elif idColumn and seqColumn:
            for x in range(2, numCells):
                id = ws.cell(row=x, column=idColumn).value
                seq = ws.cell(row=x, column=seqColumn).value
                if id and seq:
                    yield Sequence(id, seq)

        wb.close()
