from __future__ import annotations

from abc import abstractmethod
from pathlib import Path
from typing import Callable, Literal, NamedTuple

from Bio.SeqIO.FastaIO import SimpleFastaParser

from itaxotools.spart_parser import Spart as SpartParserSpart

from .handlers import FileHandler, ReadHandle, WriteHandle
from .sequences import Sequence


class Classification(NamedTuple):
    individual: str
    subset: str


class Partition(dict[str, str]):
    """Keys are individuals, values are subsets"""

    @classmethod
    def fromPath(
        cls, path: Path, handler: PartitionHandler, *args, **kwargs
    ) -> Partition:
        return handler.as_dict(path, *args, **kwargs)


class PartitionHandler(FileHandler[Classification]):
    @classmethod
    def as_dict(cls, path: Path, *args, **kwargs) -> Partition:
        spartition = Partition()
        for individual, subset in cls(path, "r", *args, **kwargs):
            spartition[individual] = subset
        return spartition

    def _open(
        self,
        path: Path,
        mode: Literal["r", "w"] = "r",
        filter: Callable[[Classification], Classification] = None,
        *args,
        **kwargs,
    ):
        self.filter = filter
        super()._open(path, mode, *args, **kwargs)

    def _iter_write(self) -> WriteHandle[Sequence]:
        raise NotImplementedError()

    def _iter_read(self, *args, **kwargs) -> ReadHandle[Classification]:
        inner_generator = self._iter_read_inner(*args, **kwargs)
        yield next(inner_generator)
        for classification in inner_generator:
            if self.filter:
                classification = self.filter(classification)
            if classification is None:
                continue
            yield classification

    @abstractmethod
    def _iter_read_inner(self, *args, **kwargs) -> ReadHandle[Classification]:
        while False:
            yield Classification()

    @staticmethod
    def subset_first_word(classification: Classification) -> Classification:
        individual, subset = classification
        try:
            first_word, rest = subset.split(" ", 1)
        except ValueError:
            print(f"Cannot split subset {subset} for individual {individual}")
            return None
        return Classification(individual, first_word)


class Tabular(PartitionHandler):
    subhandler = FileHandler.Tabular

    def _iter_read_inner(
        self,
        idHeader: str = None,
        subHeader: str = None,
        hasHeader: bool = False,
        idColumn: int = 0,
        subColumn: int = 1,
    ) -> ReadHandle[Classification]:
        if idHeader and subHeader:
            columns = (idHeader, subHeader)
            hasHeader = True
        else:
            columns = (idColumn, subColumn)

        with self.subhandler(
            self.path,
            has_headers=hasHeader,
            columns=columns,
        ) as rows:
            yield self
            for individual, subset in rows:
                # individual = sanitize(individual)
                # subset = sanitize(subset)
                yield Classification(individual, subset)


class Tabfile(Tabular, PartitionHandler):
    subhandler = FileHandler.Tabular.Tabfile


class Excel(Tabular, PartitionHandler):
    subhandler = FileHandler.Tabular.Excel


class Spart(PartitionHandler):
    def _iter_read_inner(self, spartition: str = None) -> ReadHandle[Classification]:
        spart = SpartParserSpart.fromPath(self.path)

        if spartition is None:
            spartition = spart.getSpartitions()[0]
        yield self

        for subset in spart.getSpartitionSubsets(spartition):
            for individual in spart.getSubsetIndividuals(spartition, subset):
                yield Classification(individual, subset)


class Fasta(PartitionHandler):
    def _iter_read_inner(self, separator="|") -> ReadHandle[Classification]:
        with open(self.path, "r") as handle:
            yield self
            for title, _ in SimpleFastaParser(handle):
                try:
                    individual, subset = title.split(separator, 1)
                except ValueError:
                    print(f"Could not extract partition info from fasta line: {title}")
                    continue
                yield Classification(individual, subset)

    @classmethod
    def has_subsets(self, path: Path, separator: str = "|") -> bool:
        if not separator:
            return False
        with open(path, "r") as handle:
            for title, _ in SimpleFastaParser(handle):
                data = title.split(separator, 1)
                return len(data) == 2

    @classmethod
    def guess_subset_separator(self, path: Path) -> bool | None:
        separators = "|."
        with open(path, "r") as handle:
            for title, _ in SimpleFastaParser(handle):
                for separator in separators:
                    if separator in title:
                        return separator
            return None
