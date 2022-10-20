from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterator, NamedTuple
from enum import Enum, auto
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO
from openpyxl import load_workbook

class Sequence(NamedTuple):
    id: str
    seq: str


class Sequences:
    def __init__(self, iterator: Iterator[Sequence]):
        self.iterator = iterator

    def __iter__(self):
        return self.iterator

    def __next__(self):
        return next(self.iterator)


class SequenceReader(Enum):
    FastaReader = auto()
    GenbankReader = auto()
    TabfileReader = auto()
    ExcelReader = auto()

    def __init__(self, *args):
        self._read = None

    def read(self, path: Path, **kwargs) -> Sequences:
        if self._read is None:
            raise MissingImplementation(self)
        return self._read(path, **kwargs)


class MissingImplementation(Exception):
    def __init__(self, reader: SequenceReader):
        self.reader = reader
        message = f'Missing implementation for {reader}'
        super().__init__(message)


def sequence_reader(reader: SequenceReader):
    """Decorate reader implementations"""
    def decorator(func: Callable):
        reader._read = func
    return decorator


@sequence_reader(SequenceReader.FastaReader)
def read_fasta_sequences(path: Path) -> Sequences:
    # Bio.SeqIO.FastaIO.SimpleFastaParser
    print(path)
    with open(path, 'r') as handle:
        for data in SimpleFastaParser(handle):
            yield Sequence(*data)


@sequence_reader(SequenceReader.GenbankReader)
def read_genbank_sequences(path: Path) -> Sequences:
    # Bio.GenBank.Scanner
    for data in SeqIO.parse(path, 'genbank'):
        yield Sequence(data.id, data.seq)


@sequence_reader(SequenceReader.TabfileReader)
def read_tabfile_sequences(
    path: Path,
    id: str = None,
    seq: str = None,
) -> Sequences:
    # new lazy iterator, optional headers
    # library.TabfileReader
    count = 0
    with open(path) as f:
        for line in f:
            if len(line) > 1:
                data = line.strip().split('\t', maxsplit=-1)
                print(len(data))
                if len(data) == 2:
                    id, seq = data[0], data[1]
                    if id and seq:
                        yield Sequence(id, seq)
                else:
                    try:
                        if count == 0:
                            indexId, indexSeq = data.index('seqid'), data.index('sequences')
                            count = 1
                            continue
                        yield Sequence(data[indexId], data[indexSeq])
                    except:
                        yield Sequence("0", "1")
                        # raise TypeError('Must either provide both or none of id/seq')

@sequence_reader(SequenceReader.ExcelReader)
def read_excel_sequences(
    path: Path,
    id: str = None,
    seq: str = None,
) -> Sequences:
    # new lazy iterator, optional headers
    # openpyxxl.readthedocs.io/en/stable/optimized.html
    # library.XlsxReader
    wb = load_workbook(filename=path, read_only=True)
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
