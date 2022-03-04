#!/usr/bin/env python3
from __future__ import annotations

from typing import List, Set, Any, Optional, Iterator
from abc import ABC, abstractmethod
from pathlib import Path
import re
import itertools
from enum import Enum, auto

import pandas as pd
import openpyxl

from itaxotools.DNAconvert.library.fasta import Fastafile
from itaxotools.DNAconvert.library.genbank import GenbankFile
from itaxotools.DNAconvert.library.record import Record


class InvalidPath(Exception):
    pass


class ValidFilePath(Path):
    """
    Subclass of Path that contains a path to an existing file
    """

    def __init__(self, *pathsegments: Any):
        super().__init__(*pathsegments)
        if not self.exists() or not self.is_file():
            raise InvalidPath


class _Header:

    RENAME_DICT = {
        "organism": "species",
        "specimen_identifier": "specimen_voucher",
        "specimen identifier": "specimen_voucher",
        "specimen voucher": "specimen_voucher",
    }

    def __init__(self, header: List[str]):
        self.names = list(map(_Header.rename, map(_Header.to_ascii, header)))

    @staticmethod
    def to_ascii(name: str) -> str:
        regex = r"([a-zA-Z0-9_-]+)"
        return "".join(re.findall(regex, name))

    @staticmethod
    def rename(name: str) -> str:
        if "sequence" in name:
            return "sequence"
        else:
            try:
                return _Header.RENAME_DICT[name]
            except KeyError:
                return name


class DataType(ABC):
    @abstractmethod
    @classmethod
    def from_path(cls, path: ValidFilePath, protocol: SequenceReader) -> DataType:
        pass


class DuplicatedSeqids(Exception):
    pass


class SequenceData(DataType):
    def __init__(self, path: ValidFilePath, protocol: SequenceReader):
        """
        raises `InvalidPath` if `path` does not exists
        """
        self.path = path
        self.protocol = protocol
        self.dataframe: Optional[pd.DataFrame] = None

    @classmethod
    def from_path(cls, path: ValidFilePath, protocol: SequenceReader) -> SequenceData:
        return cls(path, protocol)

    def get_dataframe(self) -> pd.DataFrame:
        if self.dataframe is None:
            self.dataframe = self.protocol.read(
                self.path, columns=["seqid", "sequence"]
            )
            try:
                self.dataframe.set_index("seqid", inplace=True, verify_integrity=True)
            except ValueError:
                raise DuplicatedSeqids()
        return self.dataframe

    def get_dataframe_chunks(self) -> Iterator[pd.DataFrame]:
        for dataframe in self.protocol.read(self.path, columns=["seqid", "sequence"]):
            try:
                dataframe.set_index("seqid", inplace=True, verify_integrity=True)
            except ValueError:
                raise DuplicatedSeqids()
            else:
                yield dataframe


class Metric(Enum):
    Uncorrelated = auto()
    Kimura2P = auto()
    JukesCantor = auto()
    UncorrelatedWithGaps = auto()

    def __str__(self) -> str:
        return {
            Metric.Uncorrelated: "pairwise uncorrected distance",
            Metric.Kimura2P: "Jukes-Cantor distance",
            Metric.JukesCantor: "Kimura-2-Parameter distance",
            Metric.UncorrelatedWithGaps: "pairwise uncorrected distance counting gaps",
        }[self]


class SequenceDistanceMatrix(DataType):
    def __init__(self, distance_matrix: pd.DataFrame):
        assert list(distance_matrix.index.names) == ["seqid1", "seqid2"]
        assert list(distance_matrix.columns) in set(Metric)
        assert distance_matrix.index.is_unique()
        self.distance_matrix = distance_matrix

    @classmethod
    def from_path(
        cls, path: ValidFilePath, protocol: SequenceReader
    ) -> SequenceDistanceMatrix:
        raise NotImplementedError


class SpeciesPartition(DataType):
    def __init__(self, species_partition: pd.DataFrame):
        assert list(species_partition.index.names) == ["seqid"]
        assert list(species_partition.columns) == ["species"]
        assert species_partition.index.is_unique()
        self.species_partition = species_partition

    @classmethod
    def from_path(
        cls, path: ValidFilePath, protocol: SequenceReader
    ) -> SpeciesPartition:
        species_partition = protocol.read(path, columns=["seqid", "species"])
        try:
            species_partition.set_index("seqid", inplace=True, verify_integrity=True)
        except ValueError:
            raise DuplicatedSeqids
        return cls(species_partition)


class SubsubspeciesPartition(DataType):
    def __init__(self, subspecies_partition: pd.DataFrame):
        assert list(subspecies_partition.index.names) == ["seqid"]
        assert list(subspecies_partition.columns) == ["subspecies"]
        assert subspecies_partition.index.is_unique()
        self.subspecies_partition = subspecies_partition

    @classmethod
    def from_path(
        cls, path: ValidFilePath, protocol: SequenceReader
    ) -> SubsubspeciesPartition:
        subspecies_partition = protocol.read(path, columns=["seqid", "subspecies"])
        try:
            subspecies_partition.set_index("seqid", inplace=True, verify_integrity=True)
        except ValueError:
            raise DuplicatedSeqids
        return cls(subspecies_partition)


class SequenceReader(ABC):
    DEFAULT_CHUNK_SIZE: int = 1000

    @abstractmethod
    @staticmethod
    def read(path: ValidFilePath, *, columns: List[str]) -> pd.DataFrame:
        pass

    @abstractmethod
    @staticmethod
    def read_chunks(
        path: ValidFilePath, *, columns: List[str], chunksize: int = DEFAULT_CHUNK_SIZE
    ) -> Iterator[pd.DataFrame]:
        pass


class ColumnsNotFound(Exception):
    def __init__(self, columns: Set[str]):
        if not isinstance(columns, set):
            raise TypeError("ColumnsNotFound argument should be a set")
        self.args = (columns,)


class TabfileReader(SequenceReader):
    @staticmethod
    def _verify_columns(path: ValidFilePath, *, columns: List[str]) -> List[str]:
        with open(path, errors="replace") as file:
            header = _Header(file.readline().rstrip("\n").split("\t"))
        if not set(columns) in set(header.names):
            raise ColumnsNotFound(set(columns) - set(header.names))
        return header.names

    @staticmethod
    def read(path: ValidFilePath, *, columns: List[str]) -> pd.DataFrame:
        names = TabfileReader._verify_columns(path, columns=columns)
        return pd.read_table(
            path,
            header=0,
            names=names,
            usecols=columns,
            na_filter=False,
            dtype=str,
            encoding_errors="replace",
        )

    @staticmethod
    def read_chunks(
        path: ValidFilePath,
        *,
        columns: List[str],
        chunksize: int = SequenceReader.DEFAULT_CHUNK_SIZE
    ) -> Iterator[pd.DataFrame]:
        names = TabfileReader._verify_columns(path, columns=columns)
        yield from pd.read_table(
            path,
            header=0,
            names=names,
            usecols=columns,
            na_filter=False,
            dtype=str,
            encoding_errors="replace",
            chunksize=chunksize,
        )


class XlsxReader(SequenceReader):
    @staticmethod
    def read(path: ValidFilePath, *, columns: List[str]) -> pd.DataFrame:
        worksheet = openpyxl.load_workbook(path).worksheets[0]
        header = _Header(list(map(lambda cell: cell.value, next(worksheet.rows))))
        if not set(columns) in set(header.names):
            raise ColumnsNotFound(set(columns) - set(header.names))
        return pd.read_excel(
            path,
            sheet_name=0,
            header=0,
            names=header.names,
            usecols=columns,
            na_filter=False,
            dtype=str,
            engine="openpyxl",
            encoding_errors="replace",
        )


class FastaReader(SequenceReader):
    @staticmethod
    def read(path: ValidFilePath, *, columns: List[str]) -> pd.DataFrame:
        if not set(columns) in {"seqid", "sequence"}:
            raise ColumnsNotFound(set(columns) - {"seqid", "sequence"})
        with open(path, errors="replace") as file:
            _, records = Fastafile.read(file)
            return pd.DataFrame(
                ([record["seqid"], record["sequence"]] for record in records()),
                columns=["seqid", "sequence"],
            )


class GenbankReader(SequenceReader):
    @staticmethod
    def select_columns(columns: List[str], record: Record) -> List[str]:
        return list(
            record[column] if column != "species" else record["organism"]
            for column in columns
        )

    @staticmethod
    def _verify_columns(path: ValidFilePath, *, columns: List[str]) -> None:
        with open(path, errors="replace") as file:
            names, _ = GenbankFile.read(file)
        if not set(columns) in set(names):
            raise ColumnsNotFound(set(columns) - set(names))

    @staticmethod
    def read(path: ValidFilePath, *, columns: List[str]) -> pd.DataFrame:
        GenbankReader._verify_columns(path, columns=columns)
        with open(path, errors="replace") as file:
            _, records = GenbankFile.read(file)
            return pd.DataFrame(
                (GenbankReader.select_columns(columns, record) for record in records()),
                columns=columns,
            )

    @staticmethod
    def read_chunks(
        path: ValidFilePath,
        *,
        columns: List[str],
        chunksize: int = SequenceReader.DEFAULT_CHUNK_SIZE
    ) -> Iterator[pd.DataFrame]:
        GenbankReader._verify_columns(path, columns=columns)
        with open(path, errors="replace") as file:
            _, records = GenbankFile.read(file)
            records = records()
            while True:
                table = pd.DataFrame(
                    (
                        GenbankReader.select_columns(columns, record)
                        for record in itertools.islice(records, chunksize)
                    ),
                    columns=columns,
                )
                if table.empty:
                    break
                yield table
