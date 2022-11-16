from __future__ import annotations

from pathlib import Path
from typing import NamedTuple
from contextlib import contextmanager

from .sequences import Sequence, Sequences
from .types import Container, Type


class Spartition(dict[str, str]):
    @classmethod
    def fromFile(cls, file: SpartitionFile, *args, **kwargs) -> Spartition:
        return file.read(*args, **kwargs)


class SpartitionFile(Type):
    """Handlers for spartition files"""

    def __init__(self, path: Path):
        self.path = path

    def read(self, *args, **kwargs) -> Spartition:
        raise NotImplementedError()


class Tabular(SpartitionFile):

    def iter_rows(self):
        raise NotImplementedError()

    @contextmanager
    def open(self):
        yield self.iter_rows()

    def read(
        self,
        idHeader: str = None,
        subsetHeader: str = None,
        hasHeader: bool = False,
        idColumn: int = 0,
        subsetColumn: int = 1,
    ) -> Spartition:

        spartition = Spartition()

        with self.open() as rows:

            # Checking headers
            if idHeader and subsetHeader:
                hasHeader = True
            if hasHeader:
                headers = next(rows)
                idColumn, subsetColumn = headers.index(idHeader), headers.index(subsetHeader)

            # Getting id and seq
            for row in rows:
                if len(row) <= 1:
                    continue
                id = row[idColumn]
                sub = row[subsetColumn]

                spartition[id] = sub
            return spartition

class Tabfile(Tabular, SpartitionFile):
    def iter_rows(self) -> iter[tuple[str, ...]]:
        with open(self.path) as file:
            for line in file:
                yield line.strip().split('\t')

class Excel(Tabular, SpartitionFile):
    def iter_rows(self) -> iter[tuple[str, ...]]:
        wb = load_workbook(filename=self.path, read_only=True)
        ws = wb.worksheets[0]
        for row in ws.iter_rows(values_only=True):
            row = list(row)
            while row and row[-1] is None:
                del row[-1]
            yield [x if x else '' for x in row]
        wb.close


class Spart(SpartitionFile):
    pass
