from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterator, NamedTuple
from enum import Enum, auto


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
    return Sequences(iter([
        Sequence('id2', 'ATG'),
        Sequence('id1', 'ATC'),
        Sequence('id3', 'ATA'),
    ]))


@sequence_reader(SequenceReader.GenbankReader)
def read_genbank_sequences(path: Path) -> Sequences:
    # Bio.GenBank.Scanner
    return Sequences(iter([]))


@sequence_reader(SequenceReader.TabfileReader)
def read_tabfile_sequences(
    path: Path,
    id: str = None,
    seq: str = None,
) -> Sequences:
    # new lazy iterator, optional headers
    # library.TabfileReader
    if id and seq:
        id_col = ...
        seq_col = ...
    elif not id and not seq:
        id_col = 0
        id_seq = 1
    else:
        raise TypeError('Must either provide both or none of id/seq')
    return Sequences(iter([]))


@sequence_reader(SequenceReader.ExcelReader)
def read_excel_sequences(
    path: Path,
    id: str = None,
    seq: str = None,
) -> Sequences:
    # new lazy iterator, optional headers
    # openpyxxl.readthedocs.io/en/stable/optimized.html
    # library.XlsxReader
    return Sequences(iter([]))
