from __future__ import annotations

from pathlib import Path
from typing import NamedTuple
from contextlib import contextmanager

from itaxotools.spart_parser import Spart as SpartParserSpart

from .sequences import Sequence, Sequences
from .types import Container, Type


class Partition(dict[str, str]):
    @classmethod
    def fromFile(cls, file: PartitionFile, *args, **kwargs) -> Partition:
        return file.get(*args, **kwargs)


class PartitionFile(Type):
    """Handlers for spartition files"""

    def __init__(self, path: Path):
        self.path = path

    def read(self, *args, **kwargs) -> iter[tuple[str, str]]:
        raise NotImplementedError()

    def get(self, *args, **kwargs) -> Partition:
        spartition = Partition()
        for individual, subset in self.read(*args, **kwargs):
            spartition[individual] = subset
        return spartition


class Tabular(PartitionFile):

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
    ) -> Partition:

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

                yield (id, sub)


class Tabfile(Tabular, PartitionFile):
    def iter_rows(self) -> iter[tuple[str, ...]]:
        with open(self.path) as file:
            for line in file:
                yield line.strip().split('\t')


class Excel(Tabular, PartitionFile):
    def iter_rows(self) -> iter[tuple[str, ...]]:
        wb = load_workbook(filename=self.path, read_only=True)
        ws = wb.worksheets[0]
        for row in ws.iter_rows(values_only=True):
            row = list(row)
            while row and row[-1] is None:
                del row[-1]
            yield [x if x else '' for x in row]
        wb.close


class Spart(PartitionFile):
    def read(self, spartition=None) -> Partition:
        spart = SpartParserSpart.fromPath(self.path)

        if spartition is None:
            spartition = spart.getSpartitions()[0]

        for subset in spart.getSpartitionSubsets(spartition):
            for individual in spart.getSubsetIndividuals(spartition, subset):
                yield (individual, subset)


class Genus(Tabfile):
    def read(self, *args, **kwargs) -> Partition:
         for id, organism in super().read(*args, **kwargs):
             genus = organism.split()[0]
             yield (id, genus)
